#!/usr/bin/env python3
"""Stage 2: harmonize sample groups and apply dataset-specific strict QC."""

from __future__ import annotations

import argparse

from uctreg.common import require_obs, set_seed, validate_unique_obs_names
from uctreg.configs import DATASET_ORDER, QC_THRESHOLDS
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root


def save_qc_plot(adata, path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = ["n_genes", "nUMIs", "mito_perc"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for axis, metric in zip(axes, metrics):
        axis.hist(adata.obs[metric], bins=50, color="skyblue", edgecolor="black")
        axis.set(title=metric, xlabel=metric, ylabel="Cell count")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def process_dataset(dataset: str, paths: ProjectPaths, overwrite: bool, plots: bool) -> None:
    source = require_input(paths.sample_merge / f"{dataset}.h5ad")
    output = prepare_output(paths.qc_filtered / f"{dataset}_qc.h5ad", overwrite)
    import scanpy as sc

    adata = sc.read_h5ad(source)
    require_obs(adata, ["id", "GSM", "n_genes", "nUMIs", "mito_perc"])
    adata.obs["group"] = adata.obs["id"].astype(str).str.extract(r"([A-Za-z]+)", expand=False)
    threshold = QC_THRESHOLDS[dataset]
    keep = (
        (adata.obs["n_genes"] > threshold["min_genes"])
        & (adata.obs["n_genes"] < threshold["max_genes"])
        & (adata.obs["nUMIs"] > threshold["min_umis"])
        & (adata.obs["nUMIs"] < threshold["max_umis"])
        & (adata.obs["mito_perc"] < threshold["max_mito"])
    )
    adata = adata[keep].copy()
    adata.obs.drop(columns=["GSM_id", "id", "sample", "source_group"], errors="ignore", inplace=True)
    validate_unique_obs_names(adata, dataset)
    if plots:
        save_qc_plot(adata, paths.qc_filtered / "figs" / f"{dataset}_qc.png")
    adata.write_h5ad(output, compression="gzip")
    print(f"[{dataset}] retained {adata.n_obs:,} cells -> {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--dataset", choices=["all", *DATASET_ORDER], default="all")
    parser.add_argument("--plots", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    paths = ProjectPaths(resolve_data_root(args.data_root))
    set_seed(args.seed)
    selected = DATASET_ORDER if args.dataset == "all" else [args.dataset]
    for dataset in selected:
        process_dataset(dataset, paths, args.overwrite, args.plots)


if __name__ == "__main__":
    main()
