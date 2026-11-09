import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# See Scanpy Tutorial
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# the total number of RNA transcript counts measured at each spatial spot
plt.hist(adata.obs["total_counts"])
plt.xlabel("Total Counts") # X-axis is total counts per spot
plt.ylabel("Number of Spots") # Y-axis is number of spots with that many counts
plt.title("Total Counts per Spot")
plt.show()
# very low counts can indicate spots with little captured RNA
# very high counts may reflect differences in tissue composition or capture efficiency

# how many different genes were detected in each spatial spot
plt.hist(adata.obs["n_genes_by_counts"], bins=40)
plt.xlabel("Number of Genes")
plt.ylabel("Number of Spots")
plt.title("Total Detected per Spot")
plt.show()
# spots with few genes detected have lower qualities

adata_filtered = adata.copy()

# To be assigned based on the threshold of the counts you want to include for further analysis
sc.pp.filter_cells(adata_filtered, min_counts=5000)
sc.pp.filter_cells(adata_filtered, max_counts=35000)


sc.pp.filter_genes(adata_filtered, min_cells=10)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

print("Original:", adata.shape)
print("Filtered:", adata_filtered.shape)

# Looking for marker associated to T-cell
print("CD3D" in adata_filtered.var_names)

sc.pl.spatial(adata_filtered, color = ["CD3D", "MS4A1", "CD68"])


# Visualization
def plot_marker_genes(adata, genes):
    """
    This function creates spatial plots to show the expression of immune marker genes across the tissue.

    Args:
        adata: The data set with the normalized gene expression and spatial coordinates.
        genes: A list of marker genes we want to visualize.
    Returns:
        A plot containing the spatial gene expression plots.
    """
    # get the x and y coordinates for each spot
    x = adata.obsm["spatial"][:, 0]
    y = adata.obsm["spatial"][:, 1]

    # create a subplot for each gene
    fig, axs = plt.subplots(1, len(genes), figsize=(5 * len(genes), 5),
                            squeeze=False)

    # loop through each gene
    for i in range(len(genes)):
        gene = genes[i]

        # get expression values for the selected gene
        expression = adata[:, gene].X.toarray().flatten()

        # plot each spot and color it based on gene expression
        scatter = axs[0, i].scatter(x, y, c=expression, cmap="viridis")

        # add title and label axes
        axs[0, i].set_title(gene)
        axs[0, i].set_xlabel("X Coordinate")
        axs[0, i].set_ylabel("Y Coordinate")

        # flip the y-axis to match tissue orientation
        axs[0, i].invert_yaxis()

        # keep the x and y axes at the same scale
        axs[0, i].set_aspect("equal")

        # add a colorbar to show expression levels
        fig.colorbar(scatter, ax=axs[0, i])

    # display all plots
    plt.show()

def plot_immune_scores(adata, cell_types):
    """
    This function creates spatial plots of the immune scores for different immune cell types across the tissue.

    Args:
        adata: The data set with the immune scores and spatial coordinates.
        score.
        cell_types: A list of immune cell scores to visualize.
    Returns:
        A spatial plot of immune scores for each cell type.
    """
    # get the x and y coordinates for each spot
    x = adata.obsm["spatial"][:, 0]
    y = adata.obsm["spatial"][:, 1]

    # create a subplot for each immune cell type
    fig, axs = plt.subplots(1, len(cell_types),
                            figsize=(5 * len(cell_types), 5),
                            squeeze=False)

    # loop through each immune cell type
    for i in range(len(cell_types)):
        cell_type = cell_types[i]

        # get immune scores for each spot
        scores = adata.obs[cell_type]

        # plot each spot and color it based on immune score
        scatter = axs[0, i].scatter(x, y, c=scores, cmap="viridis")

        # add title and label axes
        axs[0, i].set_title(cell_type)
        axs[0, i].set_xlabel("X Coordinate")
        axs[0, i].set_ylabel("Y Coordinate")

        # flip the y-axis to match tissue orientation
        axs[0, i].invert_yaxis()

        # keep the x and y axes at the same scale
        axs[0, i].set_aspect("equal")

        # add a colorbar to show immune scores
        fig.colorbar(scatter, ax=axs[0, i])

    # display all plots
    plt.show()