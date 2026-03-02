#!/usr/bin/env python
# coding: utf-8
# Quality Control for Spatial Heart Atlas

import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory to script location
os.chdir('/Users/shishirbiwalkar/Github/spatial-heart-atlas')

# Set up visual settings
esc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)
sc.settings.plot_suffix = ''

print("=" * 50)
print("SPATIAL HEART ATLAS - QUALITY CONTROL")
print("=" * 50)

#Loading Single cell data
print("\n Loading Single-Cell Data...")

adata = sc.read_10x_mtx(
    'data/scRNA_raw',  
    var_names='gene_symbols',                
    cache=True
)

#Label the data
adata.obs['sample_id'] = 'heart_1k_v3'
adata.var_names_make_unique()

print(f" Loaded {adata.n_obs} cells and {adata.n_vars} genes.")


print("\nCalculating QC metrics...")

#Calculate mitochondrial genes percentage
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.obs['pct_counts_mt'] = np.sum(adata[:, adata.var['mt']].X, axis=1) / np.sum(adata.X, axis=1) * 100

#Ribosomal genes percentage
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
adata.obs['pct_counts_ribo'] = np.sum(adata[:, adata.var['ribo']].X, axis=1) / np.sum(adata.X, axis=1) * 100

#Number of genes per cell
adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(axis=1)

#Calculate total counts per cell
adata.obs['total_counts'] = np.sum(adata.X, axis=1)

print("QC metrics calculated:")
print(f"   - Mitochondrial gene percentage (pct_counts_mt)")
print(f"   - Ribosomal gene percentage (pct_counts_ribo)")
print(f"   - Number of genes by counts (n_genes_by_counts)")
print(f"   - Total counts (total_counts)")


#Plot QC Metrics

print("\n Plotting QC metrics...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

#Total counts distribution
axes[0, 0].hist(adata.obs['total_counts'], bins=50, edgecolor='black', alpha=0.7)
axes[0, 0].set_xlabel('Total Counts')
axes[0, 0].set_ylabel('Number of Cells')
axes[0, 0].set_title('Total Counts per Cell')

#Genes by counts distribution
axes[0, 1].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7)
axes[0, 1].set_xlabel('Number of Genes')
axes[0, 1].set_ylabel('Number of Cells')
axes[0, 1].set_title('Genes Detected per Cell')

#Mitochondrial % distribution
axes[1, 0].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7, color='red')
axes[1, 0].set_xlabel('Mitochondrial Gene Percentage')
axes[1, 0].set_ylabel('Number of Cells')
axes[1, 0].set_title('Mitochondrial Gene Percentage')

#Scatter plot: Total counts vs Genes
scatter = axes[1, 1].scatter(
    adata.obs['total_counts'], 
    adata.obs['n_genes_by_counts'],
    c=adata.obs['pct_counts_mt'],
    cmap='viridis',
    alpha=0.5,
    s=5
)
axes[1, 1].set_xlabel('Total Counts')
axes[1, 1].set_ylabel('Number of Genes')
axes[1, 1].set_title('Total Counts vs Genes (colored by MT %)')
plt.colorbar(scatter, ax=axes[1, 1], label='MT %')

plt.tight_layout()
plt.savefig('figures/qc_metrics.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/qc_metrics.png")


#4. Filtering

print("\n Filtering cells and genes...")

# Store original counts
print(f"Before filtering: {adata.n_obs} cells, {adata.var.shape[0]} genes")

# Filter cells based on QC thresholds
# Keep cells with:
# - At least 200 genes detected
# - At most 25% mitochondrial reads
# - At least 1000 total counts

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs['pct_counts_mt'] < 25, :]
adata = adata[adata.obs['total_counts'] > 1000, :]

print(f"After filtering: {adata.n_obs} cells, {adata.var.shape[0]} genes")


#SAVE FILTERED DATA

print("\n Saving filtered data...")

#Save as h5ad file
adata.write('data/scRNA_qc.h5ad', compression='gzip')
print(" Saved: data/scRNA_qc.h5ad")

#Summary stats
print("\n" + "=" * 50)
print(" SUMMARY STATISTICS")
print("=" * 50)
print(f"Total cells: {adata.n_obs}")
print(f"Total genes: {adata.n_vars}")
print(f"Mean counts per cell: {adata.obs['total_counts'].mean():.2f}")
print(f"Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.2f}")
print(f"Mean mitochondrial %: {adata.obs['pct_counts_mt'].mean():.2f}%"print("=" * 50)
print("QUALITY CONTROL COMPLETE!")

