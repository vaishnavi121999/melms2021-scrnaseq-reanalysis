"""Integration + embedding logic."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import scvi

from . import utils

LOGGER = utils.setup_logger()


def integrate(config_path: str, input_path: str, output_path: str, umap_path: str) -> None:
    config = utils.load_config(config_path)
    root = Path(config["project_root"])
    input_path = utils.resolve_path(input_path, root)
    output_path = utils.resolve_path(output_path, root)
    umap_path = utils.resolve_path(umap_path, root)

    LOGGER.info("Reading QC'd object from %s", input_path)
    adata = sc.read_h5ad(input_path)

    n_top_genes = config["qc"].get("n_top_genes", 2000)
    LOGGER.info("Selecting %d highly variable genes", n_top_genes)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True, flavor="seurat_v3")

    LOGGER.info("Regressing out depth / mt / ribo and scaling")
    covariates = ["total_counts", "pct_counts_mt", "pct_counts_ribo"]
    sc.pp.regress_out(adata, covariates)
    sc.pp.scale(adata, max_value=10)

    layer = config["integration"].get("layer_counts", "counts")
    batch_key = config["integration"].get("batch_key", "Sample")
    scvi.model.SCVI.setup_anndata(adata, layer=layer, categorical_covariate_keys=[batch_key])

    latent_dim = config["integration"].get("scvi_latent_dim", 30)
    LOGGER.info("Training scVI (%d latent dims)", latent_dim)
    model = scvi.model.SCVI(adata, n_latent=latent_dim)
    model.train()
    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=1e4)

    LOGGER.info("Computing neighbors/UMAP/Leiden")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    resolution = config["integration"].get("leiden_resolution", 0.5)
    sc.tl.leiden(adata, resolution=resolution)

    fig = sc.pl.umap(adata, color=["leiden", "Sample"], wspace=0.4, frameon=False, return_fig=True)
    utils.ensure_parent(umap_path)
    fig.savefig(umap_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    utils.ensure_parent(output_path)
    adata.write_h5ad(output_path)
    LOGGER.info("Saved integrated AnnData to %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate QC'd data with scVI")
    parser.add_argument("--config", default="workflow/config.yaml")
    parser.add_argument("--input", default="results/objects/qc_filtered.h5ad")
    parser.add_argument("--output", default="results/objects/integrated.h5ad")
    parser.add_argument("--umap", default="results/figures/umap_clusters.png")
    args = parser.parse_args()
    integrate(args.config, args.input, args.output, args.umap)


if __name__ == "__main__":
    main()
