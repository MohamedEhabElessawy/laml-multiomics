"""Step 1: download the three Xena files into data/raw/.

"""

from pathlib import Path

RAW_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"

# name -> URL (fill in from the Xena dataset pages)
DATASETS: dict[str, str] = {
     "mrna": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LAML.sampleMap%2FHiSeqV2.gz",
     "methylation": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LAML.sampleMap%2FHumanMethylation450.gz",
     "survival": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLAML_survival.txt",
}


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if not DATASETS:
        print("No URLs configured. Fill in DATASETS (or place files in data/raw/ manually).")
        return
    raise NotImplementedError("Implement download once URLs are recorded.")


if __name__ == "__main__":
    main()
