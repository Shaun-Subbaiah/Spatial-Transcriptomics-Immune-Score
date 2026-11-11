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