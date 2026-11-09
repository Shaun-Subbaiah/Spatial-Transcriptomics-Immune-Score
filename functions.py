import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Vilasini Gopal Functions

def filter_data(spatial_data, min_counts=5000, max_counts=35000, min_cells=10):
    """
    This function should filter low-quality spatial spots and genes

    Parameters:
        spatial_data: contains an AnnData (Annotated Data) object that has spatial transcriptomics data
        (AnnData is a data structure in Scanpy that stores gene expression measurements and their metadata)
        min_counts: stores the minimum total counts required per spot
        max_counts: stores the maximum total coutns allowed per spot
        min_cells: stores the minimum number of spots needed to detect a gene

    Returns:
        The filtered AnnData object
    """
    # Validate the filtering parameters
    if min_counts < 0 or max_counts < 0:
        raise ValueError("Count thresholds cannot be negative")
    
    if min_counts >= max_counts:
        raise ValueError("min_counts must be smaller than max_counts")

    if min_cells < 1:
        raise ValueError("min_cells must be at least 1")
    
    # Create a copy of the original dataset
    filtered_data = spatial_data.copy()

    # Identify mitochondrial genes
    filtered_data.var["mt"] = filtered_data.var_names.str.startswith("MT-")

    # Calculate the quality control metrics for each spatial spot, including total RNA counts,
    # number of genes detected, and the percentage of mitochondrial genes
    # inplace = True stores these metrics in the AnnData object's metadata
    sc.pp.calculate_qc_metrics(filtered_data, qc_vars = ["mt"], inplace = True)

    # Filter out low-quality spatial spots based on their total RNA counts
    sc.pp.filter_cells(filtered_data, min_counts = min_counts)
    sc.pp.filter_cells(filtered_data, max_counts = max_counts)

    # Filter genes by removing genes detected in fewer than the specified number of spots
    sc.pp.filter_genes(filtered_data, min_cells = min_cells)

    return filtered_data

def normalize_data(filtered_data, target_sum = 10000):
    """
    Normalize gene expression counts across spatial spots.

    Parameters:
        filtered_data: AnnData object containing filtered spatial transcriptomics data.
        target_sum: Target total expression count for each spot after normalization.

    Returns:
        Normalized AnnData object
        (adjusts gene expression values to account for differences in total coutns across spacial spots)
    """

    if target_sum <= 0:
        raise ValueError("target_sum must be greater than zero")

    normalized_data = filtered_data.copy()

    # Normalize total counts across spacial spots
    sc.pp.normalize_total(normalized_data, target_sum = target_sum)

    # Apply a logarithmic transformation to reduce the effect of highly expressed genes
    sc.pp.log1p(normalized_data)

    return normalized_data

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