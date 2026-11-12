# LAML multi-omics: do mRNA and methylation patterns relate to AML survival?

Integrative analysis of TCGA-LAML (acute myeloid leukemia) with MOFA, following the integration workflow of ElHarouni et al., *Cancer Research* 2026;86:773-84.

## Question

Do the latent factors that MOFA learns from mRNA expression and DNA methylation relate to overall survival, independent of age?

Optional extension (not started): do the factors improve survival prediction in held-out patients compared with clinical variables alone?

## Data

UCSC Xena, TCGA LAML cohort: gene expression (HiSeq), DNA methylation (450K), curated survival, and the clinical table. Sources, file names and download links are in [`data/README.md`](data/README.md). Raw and processed data are not committed to this repository.

Cohort: patients present in all three layers with valid survival: 158 (101 deaths, 57 censored), or 146 (90 deaths) when patients with survival time 0 are dropped.

## Pipeline

Run the scripts from the project root, in order; then the notebooks from the `notebooks/` folder.

| Script | What it does | Main output |
|---|---|---|
| `scripts/01_inspect.py` | Prints sample counts, patient overlap and missing values of the raw files | printed report |
| `scripts/02_preprocess.py` | Builds the cohort and cleans mRNA and methylation | `data/processed/*.pkl`, `results/tables/cohort_summary_*.csv` |
| `scripts/03_filter_probes.py` | Removes sex-chromosome and unannotated methylation probes | `methylation_filtered_*.pkl` |
| `scripts/04_clinical.py` | Attaches age, sex, risk group, FAB subtype, white cell count; cross-checks survival times | `data/processed/clinical_*.csv` |
| `notebooks/01_explore.ipynb` | Survival curve, distributions, PCA | `results/figures/` |
| `notebooks/02_mofa.ipynb` | Trains MOFA, reads the factors (variance, clinical variables, top genes/probes), tests each factor against survival with age in the model, robustness checks | `data/processed/mofa/`, `results/tables/mofa_*.csv`, `results/figures/` |

Zero-time patients (12 deaths with survival time recorded as 0) are handled by `ZERO_TIME_RULE` in `02_preprocess.py`: `floor` keeps them at 1 day (working choice), `drop` removes them (kept for a sensitivity check).

## Environment

One environment, called `laml` (Python 3.11 recommended; MOFA needs muon, mofax, mofapy2):

```bash
conda create -n laml python=3.11
conda activate laml
pip install -r requirements.txt
```

To record exactly what is installed in your working `laml` environment: `pip freeze > requirements-freeze.txt`.

## Status

- Done: data download, cohort, cleaning, probe filtering, clinical table, exploration.
- Done (trial): first MOFA run and an unadjusted Cox model on the factors (see `docs/project_notes.md`).
- Next: run `notebooks/02_mofa.ipynb` on the real data (age-adjusted Cox per factor, factor-versus-clinical table, robustness); read the top genes and probes of the factors that matter.

Decisions, findings and open items: [`docs/project_notes.md`](docs/project_notes.md).
