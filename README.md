# Spatial Heart Atlas

Single-cell and spatial transcriptomics analysis of human heart tissue using the 10x Genomics public heart datasets. The pipeline runs quality control, preprocessing, and unsupervised clustering on both scRNA-seq and Visium spatial data, producing a cell-type atlas that can be used for downstream deconvolution and spatially-resolved expression analysis.

---

## What this project does

| Step | Script | What it produces |
|---|---|---|
| Download | `src/00_download_data.py` | Raw 10x Genomics data (`data/`) |
| Quality Control | `src/01_quality_control.py` | Filtered AnnData (`data/scRNA_qc.h5ad`), QC figures |
| Preprocessing | `src/02_preprocessing.py` | Normalized, log-transformed, scaled AnnData |
| Clustering | `src/03_clustering.py` | UMAP + Leiden clusters, cluster marker genes |

The notebook `src/01_quality_control.ipynb` mirrors the QC script for interactive exploration.

---

## Data

Public datasets from 10x Genomics, no access required:

| Dataset | Source | Size |
|---|---|---|
| Human Heart scRNA-seq (1k cells, v3 chemistry) | [10x Genomics](https://www.10xgenomics.com/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-1-standard-3-0-0) | ~filtered feature-barcode matrix |
| Human Heart Visium spatial (V1) | [10x Genomics](https://www.10xgenomics.com/datasets/human-heart-1-standard-1-1-0) | ~filtered matrix + tissue image |

Run `python src/00_download_data.py` to fetch both datasets into `data/`.

---

## Setup

```bash
# Create and activate the conda environment
conda env create -f environment.yml
conda activate spatial-heart-env

# Download raw data
python src/00_download_data.py
```

### Environment

Key dependencies (see `environment.yml` for pinned versions):

- `scanpy` — single-cell analysis
- `squidpy` — spatial transcriptomics
- `tangram-sc` — spatial deconvolution
- `jupyterlab` — interactive notebooks

---

## Running the pipeline

```bash
# Run sequentially
python src/00_download_data.py
python src/01_quality_control.py
python src/02_preprocessing.py
python src/03_clustering.py
```

Figures are written to `figures/`. Intermediate AnnData objects are written to `data/`.

---

## Outputs

| File | Description |
|---|---|
| `figures/qc_metrics.png` | Violin plots — n_genes, total counts, % mitochondrial |
| `figures/umap_clustering.png` | UMAP coloured by Leiden cluster |
| `data/scRNA_qc.h5ad` | Post-QC AnnData (cells passing filters) |
| `cache/data-scRNA_raw-matrix.h5ad` | Raw matrix cache (speeds up re-runs) |

---

## QC thresholds

Defaults used in `01_quality_control.py`:

| Filter | Value |
|---|---|
| Min genes per cell | 200 |
| Max genes per cell | 5000 |
| Max % mitochondrial reads | 20% |

Adjust these in the script for your downstream use case.

---

## Project structure

```
.
├── src/
│   ├── 00_download_data.py       Download raw 10x data
│   ├── 01_quality_control.py     QC filtering
│   ├── 01_quality_control.ipynb  Interactive QC notebook
│   ├── 02_preprocessing.py       Normalization, HVG selection, scaling
│   └── 03_clustering.py          PCA, neighbour graph, UMAP, Leiden
├── figures/                      Output plots
├── cache/                        AnnData cache files
├── environment.yml               Conda environment spec
└── .gitignore
```

---

## Author

Shishir Biwalkar
