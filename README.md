# Melms 2021 Single-Cell Lung Atlas Reanalysis

This repository turns the original exploratory notebook for the Melms et al. (Nature 2021, [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524)) COVID-19 lung atlas into a reproducible single-cell RNA-seq project. It demonstrates:

- High-dimensional omics analysis across >20 lung samples (COVID-19 vs control)
- Reusable Scanpy/scVI workflows for QC, integration, clustering, and annotation
- Sample-aware **pseudobulk differential expression** (DESeq2) plus pathway enrichment
- Snakemake automation, environment pinning, documentation, and reporting artifacts

## Repository Structure

```
melms2021-scrnaseq-reanalysis/
├── README.md
├── environment.yml              # Reproducible Conda environment (Python + R/DESeq2)
├── requirements.txt             # Pip mirror when Conda is not available
├── .gitignore                   # Excludes large data objects & intermediates
├── data/
│   ├── README.md                # Download/organization instructions
│   ├── raw/                     # (gitignored) counts from GSE171524_RAW
│   ├── processed/               # (gitignored) intermediate AnnData objects
│   └── metadata/                # Sample sheets & covariates safe for git
├── notebooks/
│   └── 01_exploration_original_notebook.ipynb
├── src/
│   ├── preprocessing.py         # QC + doublet removal + filtering
│   ├── integrate.py             # Count concatenation + scVI integration
│   ├── annotate.py              # Leiden clustering, markers, cell-type labels
│   ├── pseudobulk_de.py         # Sample-aware DE + R/DESeq2 hand-off
│   ├── enrichment.py            # Pathway enrichment utilities (gseapy)
│   └── utils.py                 # Shared helpers (logging, config, IO)
├── workflow/
│   ├── Snakefile                # End-to-end automation (Scanpy + R script)
│   └── config.yaml              # Paths, QC thresholds, comparisons of interest
├── scripts/
│   └── pseudobulk_deseq2.R      # Minimal DESeq2 runner invoked by Snakemake
├── results/
│   ├── figures/                 # Publication-style plots (UMAP, volcano, etc.)
│   ├── tables/                  # Marker lists, DE tables
│   └── objects/                 # Serialized AnnData models (gitignored)
└── docs/
    └── project_notes.md         # Running log, biological interpretation hooks
```

## Quick Start

```bash
# 1. Create the environment (Python + Scanpy + scVI + R/DESeq2)
conda env create -f environment.yml
conda activate sc-covid-atlas

# 2. Point the project to the raw GEO counts (default expects ./GSE171524_RAW)
ln -s ../GSE171524_RAW data/raw/GSE171524_RAW   # or update workflow/config.yaml

# 3. Reproduce the full workflow
snakemake --cores 4 results/tables/at2_covid_vs_control_deseq2.csv

# 4. Explore the original notebook or inspect docs/project_notes.md
jupyter lab notebooks/01_exploration_original_notebook.ipynb
```

## Workflow Highlights

| Step | Tooling | Outputs |
|------|---------|---------|
| Doublet-aware QC | Scanpy, SOLO (optional) | `results/objects/qc_filtered.h5ad` |
| Integration | scVI latent space + neighbors/UMAP | `results/objects/integrated.h5ad`, `results/figures/umap_clusters.png` |
| Annotation | Leiden clustering + curated markers | `results/tables/cell_type_markers.csv`, `results/figures/umap_celltypes.png` |
| **Pseudobulk DE (upgrade)** | Python aggregation + R/DESeq2 | `results/tables/at2_covid_vs_control_deseq2.csv`, volcano plot |
| Enrichment/reporting | gseapy enrichr + markdown notes | `results/tables/enrichment_at2.csv`, `docs/project_notes.md` |

## Result Snapshot (March 11, 2026 run)

- **Cells retained:** 108,612 cells × 2,000 HVGs after QC (top 2% gene-count trimming, mitochondrial/ribosomal filters).
- **Clusters:** 17 Leiden partitions (0–16) spanning macrophage, fibroblast, endothelial, AT1/AT2, NK/T, B/plasma, and cycling populations. Largest groups: macrophage (24k), fibroblast (20k), endothelial (15k), AT2 (15k), AT1 (7k).
- **Pseudobulk DE (AT2 COVID-19 vs control):** strong up-regulation of heat-shock/chaperone programs (HSPA1A, HSP90AA1, HSPH1, HSPD1) alongside interferon-stimulated genes; volcano plot + ranked table saved under `results/`.
- **Artifacts:** `results/figures/umap_clusters.png`, `results/figures/umap_celltypes.png`, and `results/figures/at2_covid_vs_control_volcano.png` showcase embeddings, annotations, and sample-aware DE signals suitable for MD Anderson interview packets.

---