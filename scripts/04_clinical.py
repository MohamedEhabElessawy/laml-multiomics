# Step 4: attach clinical variables (age, sex, risk group, subtype, blood counts)
# to the cohort, and cross-check the survival times against the clinical file.
#
# Run from the project root, after 02_preprocess.py has been run for each rule:
#     python scripts/04_clinical.py
#
# Output (one per zero-time rule):
#     data/processed/clinical_<rule>.csv      one row per cohort patient, same order
#                                             as the survival table (ignored by git)
#     results/tables/clinical_missing_<rule>.csv   missing values per variable
#
# Nothing is recoded here. Categories are kept exactly as written in the file,
# so we can look at them before deciding how to group them.

import pandas as pd

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
RAW = "data/raw/"
PROC = "data/processed/"
TABLES = "results/tables/"
CLINICAL_FILE = "TCGA.LAML.sampleMap_LAML_clinicalMatrix"
SURVIVAL_FILE = "survival_LAML_survival.txt"
RULES = ["drop", "floor"]

# short name -> column name in the clinical file
NUMERIC = {
    "age": "age_at_initial_pathologic_diagnosis",
    "wbc": "lab_procedure_leukocyte_result_unspecified_value",
    "bm_blast_pct": "lab_procedure_bone_marrow_blast_cell_outcome_percent_value",
    "pb_blast_pct": "lab_procedure_blast_cell_outcome_percentage_value",
}
CATEGORICAL = {
    "gender": "gender",
    "cyto_risk": "acute_myeloid_leukemia_calgb_cytogenetics_risk_category",
    "fab": "leukemia_french_american_british_morphology_code",
}

# ----------------------------------------------------------------------------
# Load the clinical table: one row per patient, indexed by the 12-character ID
# ----------------------------------------------------------------------------
clin = pd.read_csv(RAW + CLINICAL_FILE, sep="\t")
assert clin["_PATIENT"].is_unique, "some patients appear twice in the clinical file"
clin = clin.set_index("_PATIENT")
print("clinical rows:", len(clin))

# ----------------------------------------------------------------------------
# Check 1: does the clinical file agree with the survival file?
# The clinical file stores days_to_death (for patients who died) and
# days_to_last_followup (for the others). If our OS.time is right, then
#     days_to_death   if it exists,   otherwise   days_to_last_followup
# should reproduce OS.time.
# ----------------------------------------------------------------------------
surv_raw = pd.read_csv(RAW + SURVIVAL_FILE, sep="\t").set_index("_PATIENT")
both = clin.join(surv_raw[["OS", "OS.time"]], how="inner")

for col in ["days_to_death", "days_to_last_followup"]:
    both[col] = pd.to_numeric(both[col], errors="coerce")     # text like "[Not Applicable]" -> NaN
both["clinical_time"] = both["days_to_death"].fillna(both["days_to_last_followup"])

has_os = both["OS.time"].notna()
mismatch = has_os & (both["clinical_time"] != both["OS.time"])
print("\n== Check 1: clinical times vs OS.time")
print("  patients in both files:", len(both))
print("  with OS.time:", int(has_os.sum()), "| clinical time disagrees:", int(mismatch.sum()))
if mismatch.any():
    print(both.loc[mismatch, ["vital_status", "days_to_death", "days_to_last_followup",
                              "OS", "OS.time"]].head(10).to_string())

# ----------------------------------------------------------------------------
# Check 2: the patients with OS.time <= 0 (the drop-versus-floor question).
# days_to_death = 0 means the source really recorded a death at day 0.
# days_to_death empty means the 0 may be an artifact.
# ----------------------------------------------------------------------------
zero = both[both["OS.time"] <= 0]
print("\n== Check 2: patients with OS.time <= 0:", len(zero))
print(zero[["vital_status", "days_to_death", "days_to_last_followup", "OS", "OS.time"]]
      .to_string())

# ----------------------------------------------------------------------------
# Build the clinical table for each cohort, in the same patient order
# ----------------------------------------------------------------------------
for rule in RULES:
    cohort = pd.read_csv(PROC + f"survival_{rule}.csv")["PatientID"].tolist()
    absent = [p for p in cohort if p not in clin.index]
    assert not absent, f"cohort patients missing from the clinical file: {absent[:5]}"

    table = pd.DataFrame(index=pd.Index(cohort, name="PatientID"))

    for short, col in NUMERIC.items():
        raw = clin.loc[cohort, col]
        table[short] = pd.to_numeric(raw, errors="coerce").values      # non-numbers -> NaN
    for short, col in CATEGORICAL.items():
        table[short] = clin.loc[cohort, col].astype("string").str.strip().values

    missing = table.isna().sum().rename("missing")
    table.to_csv(PROC + f"clinical_{rule}.csv")
    missing.rename_axis("variable").to_csv(TABLES + f"clinical_missing_{rule}.csv")

    print(f"\n== clinical table, rule: {rule} ({len(table)} patients)")
    print("missing values per variable:")
    print(missing.to_string())

    # Print the details once, for the primary rule
    if rule == RULES[0]:
        print("\nnumeric variables:")
        print(table[list(NUMERIC)].describe().round(1).to_string())
        for short in CATEGORICAL:
            print(f"\n{short} (exact values as written in the file):")
            print(table[short].value_counts(dropna=False).to_string())
