# Spatial Transcriptomics Immune Signature Analysis

A Python based project to identify and analyze immune cell signatures in diseased tissues using spatial transcriptomics data.

## About the Project

The goal of this project is to analyze spatial transcriptomics datasets to identify immune related gene expression patterns and their spatial distribution within tissues.

We aim to develop a method for scoring immune signatures that considers both gene expression and the spatial relationships between neighboring spots.


## Dataset

**Selected dataset:** 

## Methodology

## Visualization and Summary Statistics

Four visualization functions were developed in `plotting.py`:

- `plot_marker_genes()`: This function creates spatial plots to show the expression of immune marker genes across the tissue.
- `plot_immune_scores()`: This function creates spatial plots of the immune scores for different immune cell types across the tissue.
- `compare_spatial_scores()`: This function compares the original immune scores and the immune scores that take into account neighboring spots.
- `summarize_results()`:  This function calculates summary statistics and visualizes the distribution of the original and combined immune scores for each immune cell type using histograms and boxplots.

All together, these functions allow us to visualize how incorporating neighboring spots can affect immune scores and their spatial distribution.

## Results


## AI Usage

- Shaun Subbaiah:
AI was used mainly to identify weaknesses in the functions written by Shaun. It provided a list of checks for each function and assisted with cleaning up the structure of the checks after Shaun typed them out. It also provided suggestions on improving functionality of Shaun's functions. Most significantly, the suggestion to use KDtree to find neighbors more efficiently in the find_neighbors() function, was a significant improvement to the generic looping structure. Shaun implemented the KDtree change based on guidance provided by AI. All other suggestions were minor and structural.
- Aditi Thambala:
AI was used to help construct the `plot_marker_genes()` and `plot_immune_scores()` functions. Aditi then used these functions as a foundation to develop the `compare_spatial_scores()` and `summarize_results()` functions, using AI to provide suggestions for how to code a certain part as needed.

## Contributors

- Aditi Thambala
- Shaun Subbaiah
- Vilasini Gopal
- Yufei Ma

## References

- [Initial reference paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10372073/)