import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Visualization
def plot_marker_genes(adata, genes):
    """
    This function creates spatial plots to show the expression of immune marker genes across the tissue.

    Args:
        adata: The data set with the normalized gene expression and spatial coordinates.
        genes: A list of marker genes we want to visualize.
    Returns:
        A plot containing the spatial gene expression plots.
    note: Aditi's function; used AI to generate code
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
        cell_types: A list of immune cell scores to visualize.
    Returns:
        A spatial plot of immune scores for each cell type.
    note: Aditi's function; used AI to generate code
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
        scores = adata.obs[cell_type + "_score"]

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

def compare_spatial_scores(adata, cell_types):
    """
    This function compares the original immune scores and the immune scores that take into account neighboring spots.

    Args:
        adata: The data set with the immune scores and spatial coordinates.
        cell_types: A list of immune cell types to visualize.
    Returns:
            Spatial plots comparing the original and combined immune scores for each cell type.
    note: Aditi's function; written by Aditi, but used previous plot_immune_scores() as basis
    """
    # get the x and y coordinates for each spot
    x = adata.obsm["spatial"][:, 0]
    y = adata.obsm["spatial"][:, 1]

    # create two subplots for each immune cell type
    fig, axs = plt.subplots(2, len(cell_types),
                                figsize=(5 * len(cell_types), 5),
                                squeeze=False)

    # loop through each immune cell type
    for i in range(len(cell_types)):
        cell_type = cell_types[i]

        # get the original and combined immune scores
        original = adata.obs[cell_type + "_score"]
        combined = adata.obs[cell_type + "_combined_score"]

        # get the minimum score for each type of immune score
        original_min = original.min()
        combined_min = combined.min()
        
        # get the overall minimum score
        min_score = min(original_min, combined_min)
        
        # get the maximum score for each type of immune score
        original_max = original.max()
        combined_max = combined.max()
        
        # get the overall maximum score
        max_score = max(original_max, combined_max)

        # plot the original immune scores
        scatter1 = axs[0, i].scatter(x, y, c=original,
                                    cmap="viridis",
                                    vmin=min_score,
                                    vmax=max_score)

        # plot the combined immune scores
        scatter2 = axs[1, i].scatter(x, y, c=combined,
                                    cmap="viridis",
                                    vmin=min_score,
                                    vmax=max_score)

        # add titles to both plots
        axs[0, i].set_title(cell_type + " Original Score")
        axs[1, i].set_title(cell_type + " Combined Score")

        # label axes for the original score plot
        axs[0, i].set_xlabel("X Coordinate")
        axs[0, i].set_ylabel("Y Coordinate")
        axs[0, i].invert_yaxis()
        axs[0, i].set_aspect("equal")

        # label axes for the combined score plot
        axs[1, i].set_xlabel("X Coordinate")
        axs[1, i].set_ylabel("Y Coordinate")
        axs[1, i].invert_yaxis()
        axs[1, i].set_aspect("equal")

        # add colorbars to show immune scores
        fig.colorbar(scatter1, ax=axs[0, i])
        fig.colorbar(scatter2, ax=axs[1, i])

    # display all plots
    plt.show()

def summarize_results(adata, cell_types):
    """
    This function calculates summary statistics and visualizes the distribution of the original and combined immune scores for each immune cell type.

    Args:
        adata: The data set with the original and combined immune scores.
        cell_types: A list of immune cell scores to visualize.
    Returns:
        The summary statistics for each immune cell type, and displays histograms and boxplots showing the distribution of the original and combined immune scores for each cell type.
    note: Aditi's function; written by Aditi
    """
    # create list of immune score col names
    cols = []

    # loop through each immune cell type
    for cell_type in cell_types:

        # get original and combined immune scores
        original_col = cell_type + "_score"
        combined_col = cell_type + "_combined_score"

        cols.append(original_col)
        cols.append(combined_col)

    # get the immune score
    scores = adata.obs[cols]

    # calculate summary statistics
    statistics = scores.describe()

    # create two histograms for each immune cell type
    fig, axs = plt.subplots(2, len(cell_types),
                            figsize=(5 * len(cell_types), 5),
                            squeeze=False)

    # loop through each immune cell tyope
    for i in range(len(cell_types)):
        cell_type = cell_types[i]

        # get original and combined immune scores
        original = adata.obs[cell_type + "_score"]
        combined = adata.obs[cell_type + "_combined_score"]

        # get the minimum and maximum original scores
        original_min = original.min()
        original_max = original.max()

        # get the minimum and maximum combined scores
        combined_min = combined.min()
        combined_max = combined.max()

        # get the overall minimum and maximum scores
        min_score = min(original_min, combined_min)
        max_score = max(original_max, combined_max)

        # create histograms using the same range of scores
        axs[0, i].hist(original, bins=30, range=(min_score, max_score))
        axs[1, i].hist(combined, bins=30, range=(min_score, max_score))

        # add title and label axes for original scores
        axs[0, i].set_title(cell_type + " Original Score")
        axs[0, i].set_xlabel("Immune Score")
        axs[0, i].set_ylabel("Number of Spots")

        # add title and label axes for combined scores
        axs[1, i].set_title(cell_type + " Combined Score")
        axs[1, i].set_xlabel("Immune Score")
        axs[1, i].set_ylabel("Number of Spots")

    # display the histograms
    plt.show()

    # create boxplots comparing the immune scores
    plt.figure(figsize=(5 * len(cell_types), 5))
    sns.boxplot(data=scores)

    # add title and label axes
    plt.title("Distribution of Original and Combined Immune Scores")
    plt.xlabel("Immune Score Type")
    plt.ylabel("Immune Score")

    # display the boxplot
    plt.show()

    return statistics