# Contributing

## Setup

```bash
conda env create -f environment.yml
conda activate spatial-heart-env
python src/00_download_data.py
```

## Pipeline conventions

- Each script in `src/` is numbered — run them in order (`00_` → `03_`).
- Intermediate AnnData objects are written to `data/` and read by the next step. Don't skip steps.
- Figures go to `figures/`. Don't commit large figure files unless they're the final publication-quality outputs.
- The `cache/` directory holds `.h5ad` cache files that speed up re-runs. These are git-ignored (large binaries).

## Adding a new analysis step

1. Create `src/04_your_step.py` following the existing script pattern (load from previous step's `.h5ad`, process, save).
2. If it produces figures, save to `figures/` with a descriptive name.
3. Update the pipeline table in `README.md`.

## Data

Raw data lives in `data/` (git-ignored). Never commit raw FASTQ, `.h5`, or large `.h5ad` files. The download script (`src/00_download_data.py`) is the single source of truth for how to get the data.

## Dependencies

If you add a dependency, update `environment.yml` with a pinned version. Test that `conda env create -f environment.yml` still works from scratch before opening a PR.
