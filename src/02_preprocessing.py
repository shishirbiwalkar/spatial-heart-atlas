#!/usr/bin/env python
# coding: utf-8
# Preprocessing and Scaling for Spatial Heart Atlas

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Visual settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

print("=" * 50)
print("SPATIAL HEART ATLAS - PREPROCESSING")
print("=" * 50)

# 1. Load your clean data
adata = sc.read_h5ad('data/scRNA_qc.h5ad')

print(f"Loaded Clean Data: {adata.n_obs} cells")

# 2. Normalization (Industry Standard: Log-Normalize)
# Normalize counts to 10,000 per cell (CP10k) so library size doesn't bias results
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize (log(x+1)) to handle the vast range of gene expression
sc.pp.log1p(adata)

print("✅ Normalization complete (CP10k + log1p)")

# 3. Feature Selection
# We only care about genes that change between cells (Highly Variable Genes)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

print(f"✅ Found {adata.var['highly_variable'].sum()} highly variable genes")

# Save the raw data before scaling (for differential expression later)
adata.raw = adata

# 4. Regress out effects of total counts and mitochondrial % (optional but recommended)
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

print("✅ Regressed out total counts and mitochondrial percentage")

# Scale data to unit variance (needed for PCA)
sc.pp.scale(adata, max_value=10)

print("✅ Scaling complete")

# Save the preprocessed data
adata.write('data/scRNA_preprocessed.h5ad', compression='gzip')
print("✅ Saved: data/scRNA_preprocessed.h5ad")

print("=" * 50)
print("✅ PREPROCESSING AND SCALING COMPLETE!")
print("=" * 50)

