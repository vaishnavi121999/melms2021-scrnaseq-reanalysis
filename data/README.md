# Data Management

This repository re-analyzes **Melms et al. 2021** (Nature, GSE171524 / SCP1219). Raw matrices are large and remain outside Git. Use the following layout locally:

```
data/
├── raw/
│   └── GSE171524_RAW/            # symlink or copy of GEO matrices (CSV.gz)
├── processed/
│   ├── combined.h5ad
│   ├── integrated.h5ad
│   └── annotated.h5ad
└── metadata/
    ├── sample_sheet.csv          # sample/condition mapping
    └── marker_curations.csv      # manual labels used in src/annotate.py
```

## Download instructions
1. From GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524 (Supplementary files).
2. Place or symlink the extracted folder under `data/raw/GSE171524_RAW/`.
3. Populate `metadata/sample_sheet.csv` with at least:
   - `sample_id` (e.g., `GSM5226574`)
   - `file_name` (relative to `data/raw/GSE171524_RAW`)
   - `condition` (`COVID19` vs `Control`)
   - `patient_id` (if available)
4. Update `workflow/config.yaml` if you use a different directory name.

> ⚠️ Large `.h5ad`, `.mtx`, `.loom`, `.csv.gz` files are **gitignored** to keep the repo lightweight. Document any additional data sources inside `docs/project_notes.md`.
