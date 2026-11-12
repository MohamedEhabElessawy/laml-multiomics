# Data

Raw and processed data are git-ignored. This file records where each dataset comes from so the project can be reproduced.

## Raw files (`data/raw/`)

All from UCSC Xena, TCGA LAML cohort hub. Downloaded 2026-11-06 (mRNA, methylation, survival) and later (clinical, probe annotation). The `.gz` files were decompressed and the survival file was renamed. The scripts expect these exact names.

| File in `data/raw/` | Xena dataset | Download URL | Contents |
|---|---|---|---|
| `HiSeqV2` | gene expression RNAseq (IlluminaHiSeq) | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LAML.sampleMap%2FHiSeqV2.gz | log2(x+1) RSEM; 20,530 genes x 173 samples |
| `HumanMethylation450` | Illumina Human Methylation 450K | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LAML.sampleMap%2FHumanMethylation450.gz | beta values (0-1); 485,577 probes x 194 samples; ~595 MB uncompressed |
| `survival_LAML_survival.txt` | curated survival data | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLAML_survival.txt | 200 patients; OS (event), OS.time (days) |
| `TCGA.LAML.sampleMap_LAML_clinicalMatrix` | phenotype / clinical matrix | TODO: paste the Xena URL | 200 patients; age, sex, cytogenetic risk, FAB subtype, blood counts |
| `probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy` | ID/Gene mapping for the 450K array | TODO: paste the Xena URL | 395,985 probes with chromosome and position |

## Sample barcode convention

- Patient ID: first 12 characters, e.g. `TCGA-AB-2934`.
- Sample type (characters 14-15): `01`-`09` tumor, `10`-`19` normal. All samples here are `03` (primary blood-derived tumor).

## Processed files (`data/processed/`, generated, not committed)

`<rule>` is `floor` or `drop` (how survival time 0 is handled).

| Script | Files |
|---|---|
| `02_preprocess.py` | `mrna_<rule>.pkl`, `methylation_<rule>.pkl`, `survival_<rule>.csv` (columns: `PatientID`, `Survival`, `Death`) |
| `03_filter_probes.py` | `methylation_filtered_<rule>.pkl` |
| `04_clinical.py` | `clinical_<rule>.csv` |
| `notebooks/02_mofa.ipynb` | `mofa/mofa_model.hdf5` |

## Cohort sizes

| | floor | drop |
|---|---|---|
| Patients | 158 | 146 |
| Deaths / censored | 101 / 57 | 90 / 56 |
| Genes after cleaning | 19,899 | 19,874 |
| Methylation probes after filtering | 383,572 | 383,623 |

Small result tables (`mofa_*.csv`) are committed in `results/tables/`; the trained model is not (`data/processed/` is git-ignored).
