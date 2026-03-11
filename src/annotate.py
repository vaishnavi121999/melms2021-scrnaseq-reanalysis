"""Marker detection + cell-type annotation."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from . import utils

LOGGER = utils.setup_logger()


def _export_markers(adata: sc.AnnData, path: Path) -> None:
    LOGGER.info("Ranking marker genes per Leiden cluster")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers = markers[markers["pvals_adj"] < 0.05]
    utils.ensure_parent(path)
    markers.to_csv(path, index=False)


def annotate(config_path: str, input_path: str, output_path: str, markers_path: str, umap_path: str) -> None:
    config = utils.load_config(config_path)
    root = Path(config["project_root"])
    input_path = utils.resolve_path(input_path, root)
    output_path = utils.resolve_path(output_path, root)
    markers_path = utils.resolve_path(markers_path, root)
    umap_path = utils.resolve_path(umap_path, root)

    LOGGER.info("Reading integrated object from %s", input_path)
    adata = sc.read_h5ad(input_path)

    _export_markers(adata, markers_path)

    LOGGER.info("Applying manual cell-type map")
    cell_map = config.get("cell_type_map", {})
    adata.obs["cell_type"] = adata.obs["leiden"].map(cell_map).fillna(adata.obs["leiden"])

    fig = sc.pl.umap(adata, color=["cell_type"], legend_loc="on data", frameon=False, return_fig=True)
    utils.ensure_parent(umap_path)
    fig.savefig(umap_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    utils.ensure_parent(output_path)
    adata.write_h5ad(output_path)
    LOGGER.info("Saved annotated AnnData to %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate clusters + export markers")
    parser.add_argument("--config", default="workflow/config.yaml")
    parser.add_argument("--input", default="results/objects/integrated.h5ad")
    parser.add_argument("--output", default="results/objects/annotated.h5ad")
    parser.add_argument("--markers", default="results/tables/cell_type_markers.csv")
    parser.add_argument("--umap", default="results/figures/umap_celltypes.png")
    args = parser.parse_args()
    annotate(args.config, args.input, args.output, args.markers, args.umap)


if __name__ == "__main__":
    main()
