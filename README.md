# LAML Multi-Omics: Survival Modeling with mRNA and DNA Methylation

Integrative analysis of TCGA acute myeloid leukemia (TCGA-LAML) combining gene expression (RNA-seq) and DNA methylation (450K) to model overall survival.

> Status: project skeleton. Update this section as work progresses.

## Research question

Do mRNA expression and DNA methylation profiles, alone or combined, carry prognostic information for overall survival in TCGA-LAML, beyond what either layer provides alone?

## Data (UCSC Xena, TCGA LAML hub)

| Layer | Xena dataset | Notes |
|---|---|---|
| mRNA | gene expression RNAseq, IlluminaHiSeq (n=173) | log2(x+1) RSEM normalized counts |
| Methylation | Illumina Human Methylation 450K | ~485k probes; single platform only (no 27K mixing) |
| Survival | Curated survival data | `OS` (event), `OS.time` |

Exact download links and dataset IDs are recorded in [`data/README.md`](data/README.md). Raw data is not committed to this repository.

## Design decisions

- Tumor samples only for omics layers; normals (sample code `-11`) excluded.
- Patients matched across layers by the first 12 characters of the TCGA barcode.
- All feature selection and scaling are fit inside cross-validation folds only, to avoid data leakage.
- Small cohort (~150 patients) and p >> n: regularized models and nested/repeated CV; random seeds fixed and recorded.

## Repository layout

```
data/raw/         downloaded Xena files (git-ignored)
data/processed/   cleaned, aligned matrices (git-ignored)
scripts/          numbered pipeline steps, run in order
src/laml/         reusable functions imported by scripts and notebooks
notebooks/        exploration and final report (outputs cleared before commit)
results/          figures and small tables (committed)
```

## Reproducing

```bash
conda env create -f environment.yml
conda activate laml-multiomics
python scripts/01_download.py
python scripts/02_preprocess.py
python scripts/03_intersect.py
python scripts/04_model.py
```

Scripts are stubs at this stage and are filled in step by step.

## License

To be decided. TCGA data are subject to their own terms of use.
