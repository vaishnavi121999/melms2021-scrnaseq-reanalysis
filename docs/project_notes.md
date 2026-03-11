# Project Notes

## Context
- Dataset: Melms et al. 2021 lethal COVID-19 lung atlas (GSE171524, SCP1219).
- Roles targeted: Assistant Data Scientist (Bioinformatics & Computational Biology, MD Anderson) and similar single-cell computational biology roles.
- Differentiator: reproducible Snakemake + pseudobulk DE upgrade versus a single exploratory notebook.

## Biological framing
1. Recover major epithelial, immune, and stromal populations from distal lung tissue.
2. Compare COVID-19 vs control within alveolar type 2 (AT2) and fibroblast populations using pseudobulk DESeq2.
3. Highlight pathway-level differences (e.g., interferon response, surfactant regulation).

### Latest run (2026-03-11)
- QC kept **108,612 cells** (median 2,000 HVGs) across 27 samples.
- Leiden resolution 0.5 yielded **17 clusters**; dominant labels: macrophage (24k), fibroblast (20k), endothelial (15k), AT2 (15k), AT1 (7k).
- AT2 COVID-19 vs control pseudobulk DESeq2 highlights heat-shock / proteostasis genes (HSPA1A, HSP90AA1, HSPH1, HSPD1) and interferon programs; volcano + table saved under `results/`.

## Technical framing
- QC heuristics: min 200 genes, remove top 2% n_genes outliers, <20% mitochondrial, <2% ribosomal.
- Integration: SCVI latent embedding (30 dimensions) → neighbors → UMAP/Leiden.
- Annotation: manual dictionary defined in `workflow/config.yaml:cell_type_map` plus marker search.
- Pseudobulk: counts aggregated per sample × cell type; DESeq2 run through `scripts/pseudobulk_deseq2.R` to keep statisticians happy.

## TODO backlog
- [ ] Add unit tests for `src/utils.py` and pseudobulk aggregation.
- [ ] Export HTML/Quarto summary with top plots.
- [ ] Benchmark alternative resolutions (Leiden 0.5 vs 1.0) + label transfer.
- [ ] Extend Snakemake profile for SLURM cluster execution.
