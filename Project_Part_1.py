import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# See Scanpy Tutorial
# change what we had in the functions.py into a function
# the function should take in transcriptomics datasets and return it
# user can choose what they want and which dataset to use
def get_data():
    adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
    return adata

# next functions add mitochondria data to adata.var
# calculate qc metrics data 
def calculate_qc_metrics(adata):
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    return adata

# plot the histograms from the data we have
def plot_histograms(adata):
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

# this function filter the data by selecting min and max counds and min cells
def filter_data(adata, min_counts, max_counts, min_cells):
    adata_filtered = adata.copy()
    # To be assigned based on the threshold of the counts you want to include for further analysis
    sc.pp.filter_cells(adata_filtered, min_counts=min_counts)
    sc.pp.filter_cells(adata_filtered, max_counts=max_counts)
    sc.pp.filter_genes(adata_filtered, min_cells=min_cells)    
    return adata_filtered

# normalize the data
def normalize_data(adata, target_sum):
    adata_normalized = adata.copy()
    # preserve original counts; consulted ChatGPT and it suggests to do this
    adata_normalized.layers["counts"] = (adata_normalized.X.copy()) 
    sc.pp.normalize_total(adata_normalized, target_sum=target_sum)
    sc.pp.log1p(adata_normalized)
    return adata_normalized

# combine all functions into one for convenience
def preprocess_data():
    adata = get_data()
    adata = calculate_qc_metrics(adata)
    plot_histograms(adata)
    adata = filter_data(adata)
    plot_histograms(adata) # see histograms after filtering
    adata = normalize_data(adata)
    return adata