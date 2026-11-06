# Step 0: inspect the raw Xena files before initializing anything.
#
# Run from the project root (the folder that contains README.md):
#     python scripts/00_inspect.py
#
# Nothing is saved or modified here. We only doing data inspection.

import pandas as pd


# Settings
# ----------------------------------------------------------------------------
RAW = "data/raw/"                     # folder with the three downloaded files
RUN_NA_SCAN = True                    # scanning methylation for NAs takes a few minutes
                                      # set to False to skip it

# =========================================
# A TCGA barcode looks like:  TCGA-AB-2934-03
#   first 12 characters  = the patient   (TCGA-AB-2934)
#   characters 14-15     = sample type  
# In Python, c[:12] is the patient and c[13:15] is the sample type code.
# =========================================


def describe_samples(name, columns):
    """Print how many samples/patients a file has and which sample types."""
    patients = pd.Series([c[:12] for c in columns])
    types = pd.Series([c[13:15] for c in columns])
    print(f"\n** {name} **")
    print("  samples (columns):", len(columns))
    print("  unique patients:  ", patients.nunique())
    print("  sample type codes:", types.value_counts().to_dict())
    print("  patients with more than one sample:", int(patients.duplicated().sum()))
    return set(patients)



# 1. Sample IDs of the two omics files

# In both files, rows are genes/probes and COLUMNS are samples.
# nrows=0 reads only the header line, so this is instant even for the 595 MB file.
# [1:] skips the first column name ("sample"), which is just the row-ID column.
mrna_samples = pd.read_csv(RAW + "HiSeqV2", sep="\t", nrows=0).columns[1:]
meth_samples = pd.read_csv(RAW + "HumanMethylation450", sep="\t", nrows=0).columns[1:]

mrna_patients = describe_samples("mRNA (HiSeqV2)", mrna_samples)
meth_patients = describe_samples("methylation (450K)", meth_samples)
# all patients are with a barcode (03), all cases are tumor, no control patients

# 2. Survival file (one row per patient)

surv = pd.read_csv(RAW + "survival_LAML_survival.txt", sep="\t")

print("\n== survival")
print(f"  rows: {len(surv)}", " columns:", list(surv.columns))
print("  sample type codes:", surv["sample"].str[13:15].value_counts().to_dict())
print("  unique patients:  ", surv["_PATIENT"].nunique())
print("  missing OS:", int(surv["OS"].isna().sum()), 
      " missing OS.time:", int(surv["OS.time"].isna().sum()))# OS is the event (1 = dead)

# Keep only rows that have both the event (OS) and the time (OS.time)
usable = surv.dropna(subset=["OS", "OS.time"])
print("  usable rows (OS and OS.time present):", len(usable))
print("  events (OS=1):", int((usable["OS"] == 1).sum()),
      " censored (OS=0):", int((usable["OS"] == 0).sum()))
print("  OS.time <= 0:", int((usable["OS.time"] <= 0).sum()))
print("  OS.time summary (days):")
print(usable["OS.time"].describe().round(1).to_string())

surv_patients = set(usable["_PATIENT"])


# 3. Which patients are present in ALL three files?  (set intersection: &)

cohort = mrna_patients & meth_patients & surv_patients

print("\n== patient overlap")
print("  mRNA & methylation:               ", len(mrna_patients & meth_patients))
print("  mRNA & survival:                  ", len(mrna_patients & surv_patients))
print("  methylation & survival:           ", len(meth_patients & surv_patients))
print("  all three layers (usable cohort): ", len(cohort))


# 4. Survival checks restricted to that cohort

print("\n== survival: Redaction column (all rows) ==")
# 'Redaction' flags patients that TCGA recommends excluding. Empty (NaN) = not flagged.
print(surv["Redaction"].value_counts(dropna=False).to_string())

c = surv[surv["_PATIENT"].isin(cohort)]
zero = c[c["OS.time"] <= 0]
print(f"\n=> survival QC inside the {len(cohort)}-patient cohort")
print("  events (OS=1):", int((c["OS"] == 1).sum()),
      " censored (OS=0):", int((c["OS"] == 0).sum()))
print("  OS.time <= 0:", len(zero),
      f"(events: {int((zero['OS'] == 1).sum())}, censored: {int((zero['OS'] == 0).sum())})")
print("  redacted patients in cohort:", int(c["Redaction"].notna().sum()))
print("  OS.time summary (days):")
print(c["OS.time"].describe().round(1).to_string())


# 5. mRNA values (small file, so we can read it completely)

mrna = pd.read_csv(RAW + "HiSeqV2", sep="\t", index_col=0)
print("\n== mRNA values")
print("  shape (genes x samples):", mrna.shape)
print("  NaN cells:", int(mrna.isna().sum().sum()))
print("  min:", round(float(mrna.min().min()), 2), " max:", round(float(mrna.max().max()), 2))


# 6. Methylation missing values

# The file is too big to be comfortable in memory, so we read it in chunks of
# 20,000 probes (rows). For each chunk we count how many probes contain NAs,
# then move on. Only the counts are kept.
if RUN_NA_SCAN:
    n_probes = 0
    n_with_any_na = 0
    n_with_over_20pct_na = 0
    for chunk in pd.read_csv(RAW + "HumanMethylation450", sep="\t", index_col=0,
                             chunksize=20000):
        na_fraction = chunk.isna().mean(axis=1)      # fraction of NA per probe (row)
        n_probes += len(chunk)
        n_with_any_na += int((na_fraction > 0).sum())
        n_with_over_20pct_na += int((na_fraction > 0.2).sum())
    print("\n== methylation NA scan")
    print("  probes:", n_probes, " with any NA:", n_with_any_na,
          " with >20% NA:", n_with_over_20pct_na)
else:
    print("\n(methylation NA scan skipped: RUN_NA_SCAN is False)")
