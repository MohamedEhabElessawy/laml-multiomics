# Step 2: build the analysis cohort and clean the three data layers.
#
# Run from the project root (the folder that contains README.md):
#     python scripts/02_preprocess.py
#
# Output files:
#     data/processed/mrna.pkl           genes  x patients  (ignored by git)
#     data/processed/methylation.pkl    probes x patients  (ignored by git)
#     data/processed/survival.csv       PatientID, Survival, Death
#     results/tables/cohort_summary.csv every number below (this one IS committed)

import pandas as pd

# ----------------------------------------------------------------------------
# Settings (change these here, then rerun)
# ----------------------------------------------------------------------------
RAW = "data/raw/"
OUT = "data/processed/"
TABLES = "results/tables/"

# What to do with patients whose OS.time is 0 or negative:
#   "drop"  -> remove them from the cohort
#   "floor" -> keep them and set their time to FLOOR_DAYS
ZERO_TIME_RULE = "drop"
FLOOR_DAYS = 1

CHUNKSIZE = 50000   # methylation rows read at a time (keeps memory low)

summary = {}        # we collect every number here and save it at the end

# ----------------------------------------------------------------------------
# STEP 1: survival table (one row per patient)
# ----------------------------------------------------------------------------
surv = pd.read_csv(RAW + "survival_LAML_survival.txt", sep="\t")
summary["survival_rows"] = len(surv)

# 1a. Remove patients flagged as "redacted" by TCGA (none in this file; the
#     column is completely empty, but the check makes the rule explicit).
is_redacted = surv["Redaction"].notna()
summary["redacted_dropped"] = int(is_redacted.sum())
surv = surv[~is_redacted]

# 1b. Remove patients without a survival time or event status.
n_before = len(surv)
surv = surv.dropna(subset=["OS", "OS.time"])
summary["missing_os_dropped"] = n_before - len(surv)

# 1c. Patients with OS.time <= 0 (see the discussion about drop vs floor).
is_zero = surv["OS.time"] <= 0
summary["zero_time_patients"] = int(is_zero.sum())
summary["zero_time_events"] = int((is_zero & (surv["OS"] == 1)).sum())
summary["zero_time_rule"] = ZERO_TIME_RULE

if ZERO_TIME_RULE == "drop":
    surv = surv[~is_zero]
else:  # "floor"
    surv = surv.copy()
    surv.loc[is_zero, "OS.time"] = FLOOR_DAYS

# 1d. Rename to the format we want: PatientID, Survival (days), Death (1 = died)
survival = pd.DataFrame({
    "PatientID": surv["_PATIENT"].values,
    "Survival": surv["OS.time"].astype(int).values,
    "Death": surv["OS"].astype(int).values,
})
survival = survival.drop_duplicates("PatientID")

# ----------------------------------------------------------------------------
# STEP 2: which patients are in all three files?
# ----------------------------------------------------------------------------
# Reminder: barcode characters 1-12 = patient, 14-15 = sample type.
# A sample code below "10" (01-09) is a tumor sample; 10 and above are normals.
# (Two-digit codes compare correctly as text: "03" < "10".)
mrna_header = pd.read_csv(RAW + "HiSeqV2", sep="\t", nrows=0).columns
meth_header = pd.read_csv(RAW + "HumanMethylation450", sep="\t", nrows=0).columns

mrna_tumor_cols = [c for c in mrna_header[1:] if c[13:15] < "10"]
meth_tumor_cols = [c for c in meth_header[1:] if c[13:15] < "10"]

mrna_patients = {c[:12] for c in mrna_tumor_cols}
meth_patients = {c[:12] for c in meth_tumor_cols}
surv_patients = set(survival["PatientID"])

# The cohort is the patients present in all three. Sorting gives every table
# the same patient order, which matters later when we line the layers up.
cohort = sorted(mrna_patients & meth_patients & surv_patients)
cohort_set = set(cohort)

survival = survival.set_index("PatientID").loc[cohort].reset_index()

# ----------------------------------------------------------------------------
# STEP 3: mRNA (genes x samples)
# ----------------------------------------------------------------------------
mrna = pd.read_csv(RAW + "HiSeqV2", sep="\t", index_col=0)

# keep only the tumor samples of cohort patients
wanted = [c for c in mrna_tumor_cols if c[:12] in cohort_set]
mrna = mrna[wanted]
mrna.columns = [c[:12] for c in wanted]     # rename columns to patient IDs
mrna = mrna[cohort]                         # put columns in cohort order

# Drop genes that are identical in every patient (standard deviation = 0).
# They carry no information and would break later steps such as scaling.
# This uses no survival information, so it cannot leak the outcome.
summary["genes_before"] = len(mrna)
mrna = mrna[mrna.std(axis=1) > 0]
summary["genes_after"] = len(mrna)

# ----------------------------------------------------------------------------
# STEP 4: methylation (probes x samples), read in chunks
# ----------------------------------------------------------------------------
# We load only the columns we need (cohort tumor samples), which cuts memory a
# lot, and store numbers as float32 (half the memory of the default float64).
id_col = meth_header[0]
wanted = [c for c in meth_tumor_cols if c[:12] in cohort_set]
dtypes = {c: "float32" for c in wanted}

kept_chunks = []
n_probes = 0
for chunk in pd.read_csv(RAW + "HumanMethylation450", sep="\t",
                         usecols=[id_col] + wanted, index_col=id_col,
                         dtype=dtypes, chunksize=CHUNKSIZE):
    n_probes += len(chunk)
    # dropna() removes every probe (row) that has at least one missing value.
    # In this dataset NAs belong to whole probes (the probe failed), so dropping
    # them is cleaner than guessing values.
    kept_chunks.append(chunk.dropna())

meth = pd.concat(kept_chunks)
meth.columns = [c[:12] for c in meth.columns]   # rename columns to patient IDs
meth = meth[cohort]                             # cohort order
summary["probes_before"] = n_probes
summary["probes_after"] = len(meth)

# ----------------------------------------------------------------------------
# STEP 5: safety checks, then save
# ----------------------------------------------------------------------------
# Every table must contain exactly the same patients in exactly the same order.
assert list(mrna.columns) == cohort
assert list(meth.columns) == cohort
assert list(survival["PatientID"]) == cohort

summary["cohort_patients"] = len(cohort)
summary["cohort_events"] = int(survival["Death"].sum())
summary["cohort_censored"] = int((survival["Death"] == 0).sum())

# .pkl is pandas' own fast binary format: no extra libraries needed. Only open
# .pkl files that you created yourself.
mrna.to_pickle(OUT + "mrna.pkl")
meth.to_pickle(OUT + "methylation.pkl")
survival.to_csv(OUT + "survival.csv", index=False)
pd.Series(summary, name="value").rename_axis("item").to_csv(TABLES + "cohort_summary.csv")

print("== cohort summary")
for key, value in summary.items():
    print(f"  {key}: {value}")
print("\nsaved to", OUT, "and", TABLES)
