import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from numbers import Real
from pathlib import Path
from scipy.spatial import cKDTree


# AI was used to set up this directory location info
DEFAULT_DATA_DIR = (
    Path(__file__).resolve().parent / "data" / "V1_Human_Lymph_Node"
)


def load_data(data_dir=DEFAULT_DATA_DIR):
    """Load the Visium data from a folder.

    Inputs:
    data_dir : str | pathlib.Path

    Outputs:
    anndata.AnnData
    """
    data_path = Path(data_dir).expanduser().resolve()
    if not data_path.is_dir():
        raise FileNotFoundError(f"Visium data directory does not exist: {data_path}")
    matrix_path = data_path / "filtered_feature_bc_matrix.h5"
    if not matrix_path.is_file():
        raise FileNotFoundError(f"Visium matrix file does not exist: {matrix_path}")

    adata = sc.read_visium(data_path)
    adata.var_names_make_unique()
    return adata


def calculate_qc_metrics(adata):
    """Add mitochondrial labels and QC measurements to a copy of the data.

    Inputs:
    adata : anndata.AnnData

    Outputs:
    anndata.AnnData
    """
    result = adata.copy()
    result.var["mt"] = result.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(result, qc_vars=["mt"], inplace=True)
    return result


def plot_histograms(adata, bins=40, show=True):
    """Plot histograms for total counts and detected genes per spot.

    Inputs:
    adata : anndata.AnnData
    bins : int
    show : bool

    Outputs:
    matplotlib.figure.Figure
    """
    required = ("total_counts", "n_genes_by_counts")
    missing = [column for column in required if column not in adata.obs]
    if missing:
        raise KeyError(f"QC metrics are missing from adata.obs: {missing}")
    if bins <= 0:
        raise ValueError("bins must be a positive integer")

    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(adata.obs["total_counts"].dropna(), bins=bins)
    axes[0].set(xlabel="Total Counts", ylabel="Number of Spots", title="Total Counts per Spot")
    axes[1].hist(adata.obs["n_genes_by_counts"].dropna(), bins=bins)
    axes[1].set(xlabel="Number of Genes", ylabel="Number of Spots", title="Detected Genes per Spot")
    figure.tight_layout()
    if show:
        plt.show()
    return figure


# This new version of filter_data() is a combination of versions of filter_data() from the two branches (vilasini, function/Immune_Score_Neighbors)
def filter_data(adata, min_counts=5000, max_counts=35000, min_cells=10):
    """Filter low-quality spatial spots and genes using thresholds you rpovide.
    
    Inputs:
    adata : anndata.AnnData
    min_counts : int | None
    max_counts : int | None
    min_cells : int | None

    Outputs:
    anndata.AnnData
    """

    # Validate filtering thresholds
    for name, value in {
        "min_counts": min_counts,
        "max_counts": max_counts,
        "min_cells": min_cells
    }.items():
        if value is not None and value < 0:
            raise ValueError(f"{name} cannot be negative")

    if (min_counts is not None and max_counts is not None and min_counts >= max_counts):
        raise ValueError("min_counts must be smaller than max_counts")

    if min_cells is not None and min_cells < 1:
        raise ValueError("min_cells must be at least 1")

    # Work on a copy to preserve the original dataset
    filtered_data = adata.copy()

    # Calculate QC metrics only when needed
    required_metrics = {
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt"
    }

    if not required_metrics.issubset(filtered_data.obs.columns):
        filtered_data.var["mt"] = (filtered_data.var_names.str.startswith("MT-"))

        sc.pp.calculate_qc_metrics(filtered_data, qc_vars=["mt"], percent_top=None, inplace=True)

    # Filter spatial spots
    if min_counts is not None:
        sc.pp.filter_cells(filtered_data, min_counts=min_counts)

    if max_counts is not None:
        sc.pp.filter_cells(filtered_data, max_counts=max_counts)

    # Filter genes
    if min_cells is not None:
        sc.pp.filter_genes(filtered_data, min_cells=min_cells)

    return filtered_data


def normalize_data(adata, target_sum=1e4):
    """Normalize and log-transform a copy while keeping the original counts.

    Inputs:
    adata : anndata.AnnData
    target_sum : float | None

    Outputs:
    anndata.AnnData
    """
    if target_sum is not None and target_sum <= 0:
        raise ValueError("target_sum must be positive or None")
    result = adata.copy()
    result.layers["counts"] = result.X.copy()
    sc.pp.normalize_total(result, target_sum=target_sum)
    sc.pp.log1p(result)
    return result


def plot_marker_genes(adata, genes=("CD3D", "MS4A1", "CD68"), show=True, **kwargs):
    """Show where the selected marker genes appear across the tissue.

    Inputs:
    adata : anndata.AnnData
    genes : tuple[str, ...] | list[str]
    show : bool
    kwargs : dict[str, object]

    Outputs:
    matplotlib.axes.Axes | list[matplotlib.axes.Axes] | None
    """
    genes = list(genes)
    if not genes:
        raise ValueError("At least one marker gene must be provided")
    missing = [gene for gene in genes if gene not in adata.var_names]
    if missing:
        raise KeyError(f"Marker genes not found in adata.var_names: {missing}")
    return sc.pl.spatial(adata, color=genes, show=show, **kwargs)


def run_analysis(
    data_dir=DEFAULT_DATA_DIR,
    min_counts=5000,
    max_counts=35000,
    min_cells=10,
    target_sum=1e4,
    show_plots=True,
):
    """Load, check, filter, and normalize the data in one go.

    Inputs:
    data_dir : str | pathlib.Path
    min_counts : int | None
    max_counts : int | None
    min_cells : int | None
    target_sum : float | None
    show_plots : bool

    Outputs:
    anndata.AnnData
    """
    adata = calculate_qc_metrics(load_data(data_dir))
    if show_plots:
        plot_histograms(adata)
    filtered = filter_data(adata, min_counts, max_counts, min_cells)
    if show_plots:
        plot_histograms(filtered)
    normalized = normalize_data(filtered, target_sum=target_sum)
    return normalized


def find_neighbors(adata, radius, spatial_key="spatial", include_self=False):
    """Find the spots within a set distance of each spot.

    Inputs:
    adata : anndata.AnnData
    radius : int | float
    spatial_key : str
    include_self : bool

    Outputs:
    dict[str, list[str]]
    """
    if not isinstance(radius, Real) or not np.isfinite(radius) or radius < 0:
        raise ValueError("radius must be a finite, non-negative number")
    if spatial_key not in adata.obsm:
        raise KeyError(f"Coordinates not found in adata.obsm[{spatial_key!r}]")
    if not adata.obs_names.is_unique:
        raise ValueError("adata.obs_names must be unique to identify spots")

    coordinates = np.asarray(adata.obsm[spatial_key], dtype=float)
    if (
        coordinates.ndim != 2
        or coordinates.shape[0] != adata.n_obs
        or coordinates.shape[1] < 2
    ):
        raise ValueError("Spatial coordinates must have shape (n_spots, >=2)")
    if not np.isfinite(coordinates).all():
        raise ValueError("Spatial coordinates must contain only finite values")
    if adata.n_obs == 0:
        return {}

    tree = cKDTree(coordinates)
    neighbor_indices_by_spot = tree.query_ball_point(coordinates, r=float(radius))
    spot_names = adata.obs_names
    neighbors = {}

    for spot_index, candidate_indices in enumerate(neighbor_indices_by_spot):
        # Save each distance so we can list neighbors from nearest to farthest.
        selected_neighbors = []
        for candidate_index in candidate_indices:
            if not include_self and candidate_index == spot_index:
                continue

            difference = coordinates[spot_index] - coordinates[candidate_index]
            distance = np.linalg.norm(difference)
            selected_neighbors.append((float(distance), candidate_index))

        selected_neighbors.sort()

        neighbor_names = []
        for distance, candidate_index in selected_neighbors:
            neighbor_names.append(spot_names[candidate_index])

        neighbors[spot_names[spot_index]] = neighbor_names

    return neighbors


def calculate_neighbor_scores(adata, neighbors, immune_score_key):
    """Average the immune scores of each spot's neighbors.

    Inputs:
    adata : anndata.AnnData
    neighbors : dict[str, list[str]]
    immune_score_key : str

    Outputs:
    pandas.Series
    """
    if immune_score_key not in adata.obs:
        raise KeyError(f"Immune score column {immune_score_key!r} not found in adata.obs")
    scores = pd.to_numeric(adata.obs[immune_score_key], errors="coerce")
    if not adata.obs_names.is_unique:
        raise ValueError("adata.obs_names must be unique to identify spots")
    if not isinstance(neighbors, dict):
        raise TypeError("neighbors must be the mapping returned by find_neighbors")

    unknown_spots = set(neighbors).difference(adata.obs_names)
    if unknown_spots:
        raise ValueError(
            f"Neighbor mapping contains unknown spot names: {list(unknown_spots)[:5]}"
        )

    neighbor_scores = pd.Series(
        np.nan,
        index=adata.obs_names,
        dtype=float,
        name="neighbor_score",
    )

    for spot_name in adata.obs_names:
        neighbor_names = neighbors.get(spot_name, [])
        unknown_neighbors = set(neighbor_names).difference(adata.obs_names)
        if unknown_neighbors:
            raise ValueError(
                f"Neighbors for {spot_name!r} contain unknown spots: "
                f"{list(unknown_neighbors)[:5]}"
            )

        score_total = 0.0
        score_count = 0
        for neighbor_name in neighbor_names:
            score = float(scores.loc[neighbor_name])
            if np.isfinite(score):
                score_total += score
                score_count += 1

        if score_count > 0:
            neighbor_scores.loc[spot_name] = score_total / score_count

    return neighbor_scores


def calculate_spot_scores(
    adata,
    immune_score_key,
    neighbor_score_key="neighbor_score",
    neighbor_weight=0.5,
    output_key="combined_immune_score",
):
    """Mix each spot's score with the average score of its neighbors.

    Inputs:
    adata : anndata.AnnData
    immune_score_key : str
    neighbor_score_key : str
    neighbor_weight : int | float
    output_key : str

    Outputs:
    anndata.AnnData
    """
    if immune_score_key not in adata.obs:
        raise KeyError(f"Immune score column {immune_score_key!r} not found in adata.obs")
    if neighbor_score_key not in adata.obs:
        raise KeyError(f"Neighbor score column {neighbor_score_key!r} not found in adata.obs")
    if not isinstance(neighbor_weight, Real) or not np.isfinite(neighbor_weight):
        raise ValueError("neighbor_weight must be a finite number between 0 and 1")
    if not 0 <= neighbor_weight <= 1:
        raise ValueError("neighbor_weight must be between 0 and 1 inclusive")
    if not isinstance(output_key, str) or not output_key:
        raise ValueError("output_key must be a non-empty string")

    own_scores = pd.to_numeric(adata.obs[immune_score_key], errors="coerce").astype(float)
    neighbor_scores = pd.to_numeric(
        adata.obs[neighbor_score_key], errors="coerce"
    ).astype(float)
    combined_scores = own_scores.copy()
    own_weight = 1.0 - float(neighbor_weight)

    for spot_name in adata.obs_names:
        own_score = own_scores.loc[spot_name]
        neighbor_score = neighbor_scores.loc[spot_name]

        # If there is no usable neighborhood average, preserve the own score.
        if np.isfinite(neighbor_score):
            combined_scores.loc[spot_name] = (
                own_weight * own_score
                + float(neighbor_weight) * neighbor_score
            )

    result = adata.copy()
    result.obs[output_key] = combined_scores
    return result


def create_immune_dictionary():
    """
    Creates a dictionary containing marker genes for different immune cell types.

    Returns:
        Dictionary where keys are immune cell types and values 
        are lists of their associated marker genes.
    
    This gives us a reference list of genes that can help identify 
    different immune cell populations in the lymph node tissue
    """

    immune_dictionary = {
        "T_cells": ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"],
        "B_cells": ["MS4A1", "CD79A", "CD79B", "CD19"],
        "Macrophages": ["CD68", "CD163", "LYZ", "CSF1R"],
        "Dendritic_cells": ["FCER1A", "CD1C", "CLEC10A"],
        "NK_cells": ["NKG7", "GNLY", "KLRD1", "PRF1"],
        "Neutrophils": ["FCGR3B", "CSF3R", "S100A8", "S100A9"],
        "Plasma_cells": ["SDC1", "MZB1", "JCHAIN", "XBP1"]
    }

    # checks that every immune cell type has marker genes
    for cell_type, genes in immune_dictionary.items():
        if not genes:
            raise ValueError(f"{cell_type} has no marker genes.")

    # checks for duplicate genes within each cell type
    for cell_type, genes in immune_dictionary.items():
        if len(genes) != len(set(genes)):
            raise ValueError(f"Duplicate marker genes found in {cell_type}.")

    # checks that all gene names are valid strings
    for cell_type, genes in immune_dictionary.items():
        if not all(isinstance(gene, str) and gene.strip() for gene in genes):
            raise TypeError(f"Invalid gene name found in {cell_type}.")
        
    return immune_dictionary


 def calculate_immune_scores(normalized_data, immune_dictionary):
    """

    Calculate immune cell type scores for each spatial spot using 
    the normalized expression of marker genes

    Parameters:
        normalized_data: Normalized AnnData object.
        immune_dictionary: Dictionary containing marker genes grouped by immune cell type.

    Returns:
        AnnData object with immune scores added to .obs.
        The .obs stores observations(in this case spatial spots)
    """

    if not isinstance(normalized_data, sc.AnnData):
        raise TypeError("Input must be an AnnData object.")

    scored_data = normalized_data.copy()

    for cell_type, genes in immune_dictionary.items():

        # Keep only marker genes present in the dataset
        present_genes = []

        for gene in genes:
            if gene in scored_data.var_names:
                present_genes.append(gene)

        if not present_genes:
            raise ValueError(f"No marker genes for {cell_type} found in the dataset.")

        # Extract expression values for the marker genes
        expression = scored_data[:, present_genes].X

        # Calculate average expression across marker genes for each spot
        immune_score = []

        # go through each spatial spot in the dataset
        for spot in range(expression.shape[0]):
            # start the total expression at 0
            total = 0
            # go through each marker gene for the current spot
            for gene in range(expression.shape[1]):
                # add the gene's expression value to the total
                total += expression[spot, gene]

            # divide the total by the total number of marker genes to get the average
            average = total / expression.shape[1]

            # add the average expression to the immune score list
            immune_score.append(average)

        # Store scores in the spot metadata
        scored_data.obs[f"{cell_type}_score"] = immune_score

    return scored_data

