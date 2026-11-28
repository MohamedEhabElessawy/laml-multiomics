# Project notes: decisions, findings, open items

Last updated: 2026-11-28

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
| Survival test | One Cox model per factor, alone and with age; factors in SD units; Benjamini-Hochberg across factors. Then a likelihood-ratio test of the factor added to age + cytogenetic risk, on the same patients. | About 100 deaths cannot support 15 predictors in one model; age and cytogenetic risk are the known confounders. |
| Reproducibility | After one training, `RETRAIN = False` in the notebook; all tables and figures come from that one model (seed 42) | Factor numbers and signs change between trainings, so results from different runs cannot be mixed. |
| Scope (2026-11-11) | Simplified to: do the MOFA factors relate to survival, independent of age? | The earlier plan (cross-validated prediction comparison against a clinical baseline) is kept as an optional extension. |

## Findings so far

- Cohort: 158 patients (101 deaths, 57 censored); 146 (90 / 56) without the zero-time patients.
- Median survival 577 days in the 146-patient cohort; about 30% still alive on the plateau after about 1,000 days.
- mRNA: 3,619 of 19,874 genes have a mean log2(RSEM + 1) below 1, so they are barely expressed. No patient has a globally unusual expression level.
- Methylation: beta values are bimodal as expected. Per-patient mean beta runs from 0.438 to 0.64. Five patients are above 0.56: TCGA-AB-2855 (0.640), AB-3009, AB-2966, AB-2849, AB-2877.
- PCA: no outlier patients (farthest 2.9 SD, about 3.2 expected by chance). PC1 + PC2 explain 24% of mRNA and 34% of methylation variation.
- Clinical (146 patients): age 18-88, median 57.5. Cytogenetic risk: Intermediate/Normal 88, Favorable 30, Poor 26, missing 2. FAB subtype M3 (14 patients), M1 34, M4 33, M2 32, M5 14, M0 13, M7 3, M6 2, not classified 1.

### MOFA on the real data (158 patients, seed 42, retrained model; the committed tables come from this run)

- 12 factors explain at least 2% of one layer. Total variance explained: methylation 58%, mRNA 58%.
- Factor 1 is the largest pattern (23% of methylation, 1% of mRNA); it tracks FAB subtype and age, not survival.
- **Factor 3** (shared: 5.7% methylation, 6.1% mRNA): top genes are the HOXA cluster (HOXA3, A5, A6, A7, A9, A10), HOXB6, PRDM16, NKX2-3. It tracks cytogenetic risk, FAB subtype and white cell count very strongly, and age weakly.
- **Factor 4** (shared: 5.3%, 6.3%): top genes CD34, PROM1, HMGA2, CD109, CLEC9A (immature, progenitor-like cells); low end AZU1, APOC2, IRX3, IRX5. It tracks age, cytogenetic risk, FAB subtype and white cell count.
- Survival, one factor at a time. Factor 3: HR 1.45 per SD alone (q = 0.003); with age HR 1.38 (95% CI 1.11-1.71, p = 0.003, q = 0.039); holds without floored patients (p = 0.004) and without M3 patients (p = 0.017). Factor 4: alone q = 0.065, with age p = 0.23, so it mostly reflects age and risk group. No other factor reaches q < 0.10 after age adjustment.
- **Adding cytogenetic risk** (156 patients with age and risk; 99 deaths): Factor 3 HR 1.27 (0.95-1.71), likelihood-ratio p = 0.10, q = 0.41. Its survival signal is not independent of cytogenetic risk in this cohort. The effect shrinks (HR 1.36 with age only on the same patients) but does not vanish; with 99 deaths a modest independent effect cannot be excluded.
- In this cohort, relative to Intermediate/Normal risk, Favorable has HR 0.50 (p = 0.03) and Poor HR 1.24 (p = 0.38), so the risk signal is mostly "Favorable versus the rest".
- Run to run, Factor 3's age-adjusted q moved from 0.10 (first training) to 0.04 (retrained), so borderline findings are sensitive to the random start.

## First trial notebook (2026-11-11, superseded)

A first run with 3,000 genes (scanpy) and 5,000 raw-beta probes and one Cox model containing all 15 factors suggested several survival factors. They did not hold up in the corrected analysis (no age adjustment, no multiple-testing correction, about 7 events per variable, in-sample C-index). The 355-day "short versus long term" split was dropped because it labels censored patients as short-term.

## Open items

- Separate Factor 3 from cytogenetic risk: test it inside the Intermediate/Normal group, in non-Favorable patients, and by risk group and FAB subtype.
- Check stability: retrain with two more seeds and see whether a HOX-driven factor appears each time with a similar survival effect.
- Map the top CpG sites of Factors 3 and 4 to genes with the probe annotation.
- Factor 11 / Factor 12 tracks sex (q about 0.01): look at its genes.
- SNP and cross-reactive probe filtering: not done.
- Blast percentage columns: labels unverified, not used.
- Optional: add a mutation layer (NPM1, FLT3, etc.); validate in an external AML cohort (availability not checked).

## Stated from background knowledge, not verified

- FAB M3 corresponds to acute promyelocytic leukemia.
- HOXA-cluster expression is high in NPM1-mutated and KMT2A-rearranged AML and low in favorable-risk subtypes (promyelocytic, core-binding-factor); high HOXA9 is generally linked to a poorer outcome. No mutation data here to check this.
- CD34, PROM1 and HMGA2 high expression marks immature, stem-like leukemia cells.
- IDH1/IDH2-mutant AML tends to show hypermethylation.
- `OS.time` is in days.
- The Xena download URLs still resolve.
