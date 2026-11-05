"""Step 4: survival modeling.

Planned design:
  - baselines: mRNA only, methylation only, concatenated (early integration)
  - regularized Cox models (elastic net) and/or random survival forest
  - feature selection and scaling fit inside CV folds only (no leakage)
  - repeated stratified CV on the event indicator; metric: concordance index
  - fixed random seed recorded in results/tables
"""

SEED = 42


def main() -> None:
    raise NotImplementedError("Implement after 03_intersect.py.")


if __name__ == "__main__":
    main()
