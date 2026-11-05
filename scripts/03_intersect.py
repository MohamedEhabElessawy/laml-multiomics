"""Step 3: align the three layers on shared patients.

Planned steps:
  - collapse each omics layer to one column per patient (src.laml.barcodes)
  - intersect patient IDs across mRNA, methylation and survival
  - drop patients with missing survival time or event
  - save aligned matrices with identical patient order
  - write the final patient count to results/tables/cohort_summary.csv
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROC_DIR = ROOT / "data" / "processed"
TABLES_DIR = ROOT / "results" / "tables"


def main() -> None:
    raise NotImplementedError("Implement after 02_preprocess.py.")


if __name__ == "__main__":
    main()
