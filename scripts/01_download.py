"""Step 1: download the three Xena files into data/raw/.

TODO: paste the exact Xena download URLs into data/README.md, then list them
in DATASETS below. Files are saved with their original names and skipped if
they already exist.
"""

from pathlib import Path

RAW_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"

# name -> URL (fill in from the Xena dataset pages)
DATASETS: dict[str, str] = {
    # "mrna": "https://.../HiSeqV2.gz",
    # "methylation": "https://.../HumanMethylation450.gz",
    # "survival": "https://.../LAML_survival.txt",
}


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if not DATASETS:
        print("No URLs configured. Fill in DATASETS (or place files in data/raw/ manually).")
        return
    raise NotImplementedError("Implement download once URLs are recorded.")


if __name__ == "__main__":
    main()
