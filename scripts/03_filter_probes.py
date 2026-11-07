# Step 3: remove unreliable methylation probes.
#
# Run from the project root, after 02_preprocess.py has been run once for each
# zero-time rule:
#     python scripts/03_filter_probes.py
#
# What this script removes (none of it uses survival, so it cannot leak the outcome):
#   1. probes that are not in the annotation file (we cannot check where they are)
#   2. probes on the sex chromosomes X and Y (they separate men from women,
#      not leukemia biology)
#
# NOT done here: SNP-affected and cross-reactive probes. They need separate
# published exclusion lists and will be a later step.
#
# Output (one pair per zero-time rule, ignored by git except the summary table):
#     data/processed/methylation_filtered_<rule>.pkl
#     results/tables/probe_filter_summary_<rule>.csv

import pandas as pd


# Settings
# ----------------
RAW = "data/raw/"
PROC = "data/processed/"
TABLES = "results/tables/"
ANNOTATION = "probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy"

RULES = ["drop", "floor"]       # remove one if you have not created its files yet
SEX_CHROMS = ["X", "Y"]

# ----------------------------------------------------------------------------
# Load the annotation: one row per probe, with its chromosome and position.
# The first column is called "#id" in the file, so we use it as the row index.
# ----------------------------------------------------------------------------
annot = pd.read_csv(RAW + ANNOTATION, sep="\t", index_col=0)
annot.index.name = "probe"
assert annot.index.is_unique, "annotation has duplicated probe IDs"

# Chromosome labels look like "chr7"; remove the "chr" so they read "7", "X", "Y".
annot["chrom"] = annot["chrom"].astype(str).str.replace("chr", "", regex=False)

print("annotation probes:", len(annot))
print("chromosome labels in the annotation file:")
print(annot["chrom"].value_counts().to_string())


# Filter the methylation matrix for each zero-time rule

for rule in RULES:
    meth = pd.read_pickle(PROC + "methylation_" + rule + ".pkl")   # probes x patients
    summary = {"zero_time_rule": rule, "probes_start": len(meth)}

    # 1. Keep only probes that appear in the annotation.
    #    isin() gives True/False per probe; ~ means "not".
    mapped = annot.index[annot["chrom"] != "nan"] 
    in_annotation = meth.index.isin(mapped)
    summary["probes_not_in_annotation"] = int((~in_annotation).sum())
    meth = meth[in_annotation]

    # 2. Look up the chromosome of each remaining probe and drop X and Y.
    chrom = annot.loc[meth.index, "chrom"]
    on_sex_chrom = chrom.isin(SEX_CHROMS).values
    summary["probes_on_sex_chromosomes"] = int(on_sex_chrom.sum())
    meth = meth[~on_sex_chrom]

    summary["probes_final"] = len(meth)

    # Information only: which chromosome labels are left? (Anything that is not
    # 1-22, such as "." or "M", is unusual and worth a look.)
    left = annot.loc[meth.index, "chrom"].value_counts()
    odd = left[~left.index.isin([str(i) for i in range(1, 23)])]
    summary["probes_with_unusual_chrom_label"] = int(odd.sum())

    meth.to_pickle(PROC + "methylation_filtered_" + rule + ".pkl")
    pd.Series(summary, name="value").rename_axis("item").to_csv(
        TABLES + "probe_filter_summary_" + rule + ".csv")

    print(f"\n== rule: {rule}")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    if len(odd) > 0:
        print("  unusual chromosome labels still present:", odd.to_dict())
