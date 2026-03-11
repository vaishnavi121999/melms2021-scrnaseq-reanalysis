"""Lightweight wrapper around gseapy for downstream interpretation."""
from __future__ import annotations

import argparse
from pathlib import Path

import gseapy as gp
import pandas as pd

from . import utils

LOGGER = utils.setup_logger()


def enrich(results_csv: str, direction: str, library: str, top_n: int, output_csv: str) -> None:
    df = pd.read_csv(results_csv)
    column = "log2FoldChange" if "log2FoldChange" in df.columns else "lfc"
    ascending = direction == "down"
    genes = df.sort_values(column, ascending=ascending).head(top_n)[\"gene\"].tolist()
    LOGGER.info("Running enrichr on %d genes (%s-regulated)", len(genes), direction)
    enr = gp.enrichr(gene_list=genes, gene_sets=[library])
    res = enr.results
    utils.ensure_parent(output_csv)
    res.to_csv(output_csv, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run enrichr on DE results")
    parser.add_argument("results_csv")
    parser.add_argument("--direction", choices=["up", "down"], default="up")
    parser.add_argument("--library", default="GO_Biological_Process_2021")
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--output", default="results/tables/enrichment.csv")
    args = parser.parse_args()
    enrich(args.results_csv, args.direction, args.library, args.top_n, args.output)


if __name__ == "__main__":
    main()
