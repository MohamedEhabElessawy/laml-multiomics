"""Step 2: load each layer, keep tumor samples, basic cleaning.

Planned steps (fill in after inspecting the raw files in a notebook):
  mRNA:        load matrix, orient as genes x samples, keep tumor columns.
  Methylation: load matrix, keep tumor columns, drop probes with many NAs,
               drop probes on sex chromosomes. SNP / cross-reactive probe
               filtering comes later, once the probe annotation is chosen.
  Survival:    keep PatientID, OS.time, OS; rename to Survival / Death.
Save outputs to data/processed/.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"


def main() -> None:
    PROC_DIR.mkdir(parents=True, exist_ok=True)
    raise NotImplementedError("Implement after inspecting raw files.")


if __name__ == "__main__":
    main()
