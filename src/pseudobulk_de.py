"""Sample-aware pseudobulk differential expression helper."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc

from . import utils

LOGGER = utils.setup_logger()


def _resolve_layer(adata: sc.AnnData) -> np.ndarray:
    if "counts" in adata.layers:
        return adata.layers["counts"]
    if adata.raw is not None:
        return adata.raw.X
    return adata.X


def _aggregate_counts(adata: sc.AnnData, cell_type: str, min_cells: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    subset = adata[adata.obs["cell_type"] == cell_type].copy()
    if subset.n_obs < min_cells:
        raise ValueError(f"Not enough {cell_type} cells ({subset.n_obs}) for pseudobulk")
    counts = _resolve_layer(subset)
    samples = subset.obs[["Sample", "condition"]].drop_duplicates()
    matrices = []
    metadata_rows = []
    for _, row in samples.iterrows():
        sample_id = row["Sample"]
        mask = (subset.obs["Sample"] == sample_id).values
        sample_counts = counts[mask]
        summed = np.asarray(sample_counts.sum(axis=0)).ravel()
        matrices.append(summed)
        metadata_rows.append({"sample_id": sample_id, "condition": row["condition"]})
    counts_df = pd.DataFrame(matrices, columns=subset.var_names, index=[m["sample_id"] for m in metadata_rows])
    metadata_df = pd.DataFrame(metadata_rows)
    return counts_df, metadata_df


def run_deseq2(config_path: str, input_path: str, comparison: str, counts_out: str, metadata_out: str,
               results_out: str, volcano_out: str) -> None:
    config = utils.load_config(config_path)
    root = Path(config["project_root"])
    input_path = utils.resolve_path(input_path, root)
    counts_out = utils.resolve_path(counts_out, root)
    metadata_out = utils.resolve_path(metadata_out, root)
    results_out = utils.resolve_path(results_out, root)
    volcano_out = utils.resolve_path(volcano_out, root)

    comp = next((c for c in config.get("comparisons", []) if c["name"] == comparison), None)
    if comp is None:
        raise ValueError(f"Comparison '{comparison}' not found in config")

    adata = sc.read_h5ad(input_path)
    cell_type = comp["cell_type"]
    min_cells = comp.get("min_cells", 100)

    LOGGER.info("Aggregating pseudobulk counts for %s", cell_type)
    counts_df, metadata_df = _aggregate_counts(adata, cell_type, min_cells)

    utils.ensure_parent(counts_out)
    counts_df.to_csv(counts_out)
    metadata_df.to_csv(metadata_out, index=False)
    counts_per_group = metadata_df["condition"].value_counts()
    for label, count in counts_per_group.items():
        if count < 2:
            LOGGER.warning("Only %d replicates for condition %s", count, label)

    r_script = utils.resolve_path("scripts/pseudobulk_deseq2.R", root)
    cmd = [
        "Rscript",
        str(r_script),
        "--counts", str(counts_out),
        "--metadata", str(metadata_out),
        "--case", comp["case"],
        "--control", comp["control"],
        "--output", str(results_out),
        "--volcano", str(volcano_out),
    ]
    LOGGER.info("Running DESeq2 via %s", " ".join(str(x) for x in cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pseudobulk DE via DESeq2")
    parser.add_argument("--config", default="workflow/config.yaml")
    parser.add_argument("--input", default="results/objects/annotated.h5ad")
    parser.add_argument("--comparison", default="at2_covid_vs_control")
    parser.add_argument("--counts", default="results/tables/at2_counts.csv")
    parser.add_argument("--metadata", default="results/tables/at2_metadata.csv")
    parser.add_argument("--output", default="results/tables/at2_covid_vs_control_deseq2.csv")
    parser.add_argument("--volcano", default="results/figures/at2_covid_vs_control_volcano.png")
    args = parser.parse_args()
    run_deseq2(args.config, args.input, args.comparison, args.counts, args.metadata, args.output, args.volcano)


if __name__ == "__main__":
    main()
