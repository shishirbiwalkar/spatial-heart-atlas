#!/usr/bin/env python
# coding: utf-8
# Dimensionality Reduction and Clustering for Spatial Heart Atlas

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os

# Set working directory to script location
os.chdir('/Users/shishirbiwalkar/Github/spatial-heart-atlas')

# Visual settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

print("=" * 50)
print("SPATIAL HEART ATLAS - DIMENSIONALITY REDUCTION & CLUSTERING")
print("=" * 50)

# Load preprocessed data
adata = sc.read_h5ad('data/scRNA_preprocessed.h5ad')
print(f"Loaded Preprocessed Data: {adata.n_obs} cells, {adata.n_vars} genes")

# 1. Principal Component Analysis (PCA)
# Reduces data from 20,000 genes to 50 essential components
sc.tl.pca(adata, svd_solver='arpack')
print("✅ PCA complete (50 components)")

# 2. Neighborhood Graph
# Who is similar to whom? (k-Nearest Neighbors)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
print("✅ Neighborhood graph complete")

# 3. UMAP (Visualization)
# Projects the graph into 2D so we can see it
sc.tl.umap(adata)
print("✅ UMAP complete")

# 4. Clustering (Leiden Algorithm)
# Groups cells into numbered clusters (0, 1, 2, etc.)
# resolution=0.5 gives fewer, broader clusters. 1.0 gives more granular ones.
sc.tl.leiden(adata, resolution=0.5)
print("✅ Leiden clustering complete")

# Plot the result
sc.pl.umap(adata, color=['leiden'], title="Unannotated Clusters", legend_loc='on data', 
           save='_clustering.png')
print("✅ Saved: figures/umap_clustering.png")

# Save the clustered data
adata.write('data/scRNA_clustered.h5ad', compression='gzip')
print("✅ Saved: data/scRNA_clustered.h5ad")

# Print cluster summary
print("\n" + "=" * 50)
print("CLUSTER SUMMARY")
print("=" * 50)
print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
print(f"Cluster distribution:")
print(adata.obs['leiden'].value_counts().sort_index())
print("=" * 50)
print("✅ DIMENSIONALITY REDUCTION & CLUSTERING COMPLETE!")
print("=" * 50)

