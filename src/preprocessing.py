"""QC + preprocessing for Melms 2021 single-cell counts."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from . import utils

LOGGER = utils.setup_logger()


def _load_counts(sample_row: pd.Series, raw_dir: Path) -> sc.AnnData:
    file_name = sample_row["file_name"]
    path = raw_dir / file_name
    if not path.exists():
        raise FileNotFoundError(f"Missing matrix: {path}")
    LOGGER.info("Reading %s", path)
    adata = sc.read_csv(path).T
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    adata.obs["Sample"] = sample_row["sample_id"]
    adata.obs["condition"] = sample_row.get("condition", "unknown")
    if "patient_id" in sample_row:
        adata.obs["patient_id"] = sample_row["patient_id"]
    return adata


def _mark_qc_features(adata: sc.AnnData) -> None:
    genes = adata.var_names.str.upper()
    adata.var["mt"] = genes.str.startswith("MT-") | genes.str.startswith("MT.")
    adata.var["ribo"] = genes.str.startswith("RPL") | genes.str.startswith("RPS")


def _filter_cells(adata: sc.AnnData, config: dict) -> sc.AnnData:
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)
    min_genes = config["qc"].get("min_genes", 200)
    pct_mt_max = config["qc"].get("pct_mt_max", 20)
    pct_ribo_max = config["qc"].get("pct_ribo_max", 2)
    upper_quant = config["qc"].get("upper_gene_quantile", 0.98)

    LOGGER.info("Applying filters: min_genes=%s pct_mt<%s pct_ribo<%s quantile<%s",
                min_genes, pct_mt_max, pct_ribo_max, upper_quant)
    upper_lim = np.quantile(adata.obs["n_genes_by_counts"], upper_quant)
    mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["pct_counts_mt"] < pct_mt_max)
        & (adata.obs["pct_counts_ribo"] < pct_ribo_max)
        & (adata.obs["n_genes_by_counts"] < upper_lim)
    )
    LOGGER.info("Keeping %.2f%% of cells after QC", mask.mean() * 100)
    return adata[mask].copy()


def _filter_genes(adata: sc.AnnData, config: dict) -> sc.AnnData:
    min_cells = config["qc"].get("min_cells", 3)
    LOGGER.info("Filtering genes detected in < %s cells", min_cells)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata


def _optional_doublet_removal(adata: sc.AnnData, config: dict) -> sc.AnnData:
    doublet_cfg = config.get("doublet_detection", {})
    if not doublet_cfg.get("enabled", False):
        return adata
    try:
        import scvi
    except ImportError:  # pragma: no cover - optional dependency
        LOGGER.warning("scvi-tools not installed; skipping doublet detection")
        return adata
    LOGGER.info("Running SOLO doublet detection")
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train(max_epochs=doublet_cfg.get("max_epochs", 200))
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    threshold = doublet_cfg.get("decision_threshold", 0.5)
    mask = df["doublet"] < threshold
    LOGGER.info("Removing %d predicted doublets", (~mask).sum())
    return adata[mask.to_numpy()].copy()


def preprocess(config_path: str, output_path: str) -> None:
    config = utils.load_config(config_path)
    root = Path(config["project_root"])
    raw_dir = root / config["raw_data_dir"]
    sample_sheet = root / config["sample_sheet"]

    samples = utils.read_sample_sheet(sample_sheet)
    adatas: List[sc.AnnData] = []
    for _, row in samples.iterrows():
        adata = _load_counts(row, raw_dir)
        adatas.append(adata)
    LOGGER.info("Concatenating %d samples", len(adatas))
    combined = sc.concat(adatas, join="outer")
    if not sparse.issparse(combined.X):
        combined.X = sparse.csr_matrix(combined.X)
    combined.layers["counts"] = combined.X.copy()

    _mark_qc_features(combined)
    combined = _filter_genes(combined, config)
    combined = _filter_cells(combined, config)

    LOGGER.info("Normalizing and log-transforming")
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    combined.raw = combined

    combined = _optional_doublet_removal(combined, config)

    utils.ensure_parent(output_path)
    combined.write_h5ad(output_path)
    LOGGER.info("Saved QC'd AnnData to %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="QC + preprocess Melms 2021 counts")
    parser.add_argument("--config", default="workflow/config.yaml")
    parser.add_argument("--output", default="results/objects/qc_filtered.h5ad")
    args = parser.parse_args()
    preprocess(args.config, args.output)


if __name__ == "__main__":
    main()
