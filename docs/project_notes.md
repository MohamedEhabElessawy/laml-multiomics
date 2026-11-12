# Project notes: decisions, findings, open items

Last updated: 2026-11-11

## Decisions

| Decision | Choice | Reason |
|---|---|---|
| Data source | UCSC Xena (HiSeqV2, 450K, curated survival, clinical matrix) | Pre-merged matrices; GDC would need hundreds of per-sample files. Xena legacy data is not GDC-harmonized. |
| Cohort | Patients in all three layers with valid survival; tumor samples only | A model needs the same patients in every layer. |
| Survival time 0 | Working choice: `floor` (keep at 1 day). `drop` kept for a sensitivity check. | The clinical file records `days_to_death = 0` explicitly for 12 of the 13 patients, so the deaths are certain and only the timing is doubtful. The clinical file also reproduces `OS.time` for all 186 patients. |
| Methylation probes | Removed: probes with missing values (92,425 across all samples), unannotated (68), unknown chromosome (35), chrX/chrY (9,691) | Failed probes and sex-chromosome probes are not reliable signal for this question. SNP-affected and cross-reactive probes are not yet removed. |
| Integration method | MOFA (muon / mofapy2) instead of PCA, following ElHarouni et al. | MOFA shows which patterns are shared between layers and which belong to one layer. |
| MOFA input | 5,000 most variable genes and 5,000 most variable probes (variance taken on M-values), `scale_views=True`, 15 starting factors, Gaussian likelihoods, one run (seed 42); a factor is kept if it explains at least 2% of one layer | Equal feature counts and scaled views so neither layer dominates; settings follow the paper where muon allows. |
| MOFA code | Lives in `notebooks/02_mofa.ipynb` (replaces `scripts/05_run_mofa.py`) | The trial notebook already ran on the real data; one implementation, not two. |
| Survival test | One Cox model per factor, alone and with age; factors in SD units; Benjamini-Hochberg across factors | About 100 deaths cannot support 15 predictors in one model; age is the main known confounder. |
| Scope (2026-11-11) | Simplified to: do the MOFA factors relate to survival, independent of age? | The earlier plan (cross-validated prediction comparison against a clinical baseline) is kept as an optional extension. |

## Findings so far

- Cohort: 158 patients (101 deaths, 57 censored); 146 (90 / 56) without the zero-time patients.
- Median survival 577 days in the 146-patient cohort; about 30% still alive on the plateau after about 1,000 days.
- mRNA: 3,619 of 19,874 genes have a mean log2(RSEM + 1) below 1, so they are barely expressed. No patient has a globally unusual expression level.
- Methylation: beta values are bimodal as expected. Per-patient mean beta runs from 0.438 to 0.64. Five patients are above 0.56: TCGA-AB-2855 (0.640), AB-3009, AB-2966, AB-2849, AB-2877.
- PCA: no outlier patients (farthest 2.9 SD, about 3.2 expected by chance). PC1 + PC2 explain 24% of mRNA and 34% of methylation variation.
- Clinical (146 patients): age 18-88, median 57.5. Cytogenetic risk: Intermediate/Normal 88, Favorable 30, Poor 26, missing 2. FAB subtype M3 (14 patients), M1 34, M4 33, M2 32, M5 14, M0 13, M7 3, M6 2, not classified 1.
- MOFA, first real run (158 patients): two random starts gave the same ELBO and 12 factors above 2%, so the solution is stable.

## Trial notebook (2026-11-11, `trial_mofa.ipynb`, superseded by `02_mofa.ipynb`)

First end-to-end MOFA run on 158 patients: 3,000 genes (scanpy highly-variable-genes) and 5,000 probes (beta values, unfiltered methylation file), 15 factors, converged.

- Variance explained: mRNA 57%, methylation 59% in total. Factors 1, 2 and 8 are mostly methylation; 4 and 6 mostly mRNA; 3 and 5 are shared.
- Factor 3 (mRNA) is led by HOX genes (HOXA10, HOXB4, HOXB5, HOXB9, PBX3), a known AML expression pattern.
- Cox model with all 15 factors together: Factors 3, 6, 7, 14 and 15 had p < 0.05; in-sample C-index 0.70. Not reliable as it stands: no age adjustment, no multiple-testing correction, about 7 events per variable, factors on different scales, and the C-index is measured on the training patients. Factor 15 explains only 0.02% of methylation and 1.6% of mRNA yet was the strongest hit, which is a warning sign.
- The "short-term versus long-term" split at 355 days labels censored patients as short-term and is not used any more.

## Open items

- Run `notebooks/02_mofa.ipynb` on the real data and read: factor-versus-clinical table, age-adjusted Cox per factor, robustness table, top-weighted genes and probes. Check whether Factor 3 (HOX) still stands after age adjustment and whether any factor is just sex, age or FAB.
- Check the results with and without the zero-time patients, the five hypermethylated patients, and the M3 patients.
- SNP and cross-reactive probe filtering: not done.
- Blast percentage columns: labels unverified, not used.
- Optional: add a mutation layer; validate in an external AML cohort (availability not checked).

## Stated from background knowledge, not verified

- FAB M3 corresponds to acute promyelocytic leukemia.
- IDH1/IDH2-mutant AML tends to show hypermethylation.
- `OS.time` is in days.
- The Xena download URLs still resolve.
