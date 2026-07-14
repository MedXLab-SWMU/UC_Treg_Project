#!/usr/bin/env python3
"""Stage 3: merge the 11 colon-tissue datasets and create counts/raw layers."""

from __future__ import annotations

import argparse
import math
import re

from uctreg.common import set_seed, validate_unique_obs_names
from uctreg.configs import DATASET_ORDER
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root


def minimum_dataset_count(n_datasets: int, fraction: float) -> int:
    if n_datasets < 1:
        raise ValueError("n_datasets must be positive")
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    return max(1, math.ceil(n_datasets * fraction))


def dataset_id(filename: str) -> str:
    match = re.search(r"(GSE\d+|SCP\d+)", filename, flags=re.IGNORECASE)
    return match.group(1).upper() if match else filename.removesuffix("_qc.h5ad")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--min-fraction", type=float, default=0.8)
    args = parser.parse_args()
    set_seed(args.seed)

    paths = ProjectPaths(resolve_data_root(args.data_root))
    tissue_datasets = [dataset for dataset in DATASET_ORDER if "pbmc" not in dataset.lower()]
    files = [paths.qc_filtered / f"{dataset}_qc.h5ad" for dataset in tissue_datasets]
    for path in files:
        require_input(path, "tissue QC dataset")

    import numpy as np
    import scanpy as sc
    from scipy.sparse import issparse

    adatas = []
    for path in files:
        adata = sc.read_h5ad(path)
        validate_unique_obs_names(adata, path.name)
        adata.obs["GSE"] = dataset_id(path.name)
        adata.obs["GSE"] = adata.obs["GSE"].astype("category")
        adatas.append(adata)
        print(f"Loaded {path.name}: {adata.n_obs:,} cells")

    gene_sets = [set(adata.var_names) for adata in adatas]
    all_genes = set().union(*gene_sets)
    required = minimum_dataset_count(len(adatas), args.min_fraction)
    common_genes = sorted(
        gene for gene in all_genes if sum(gene in genes for genes in gene_sets) >= required
    )
    if not common_genes:
        raise ValueError("No genes passed the cross-dataset occurrence threshold")

    merged = sc.concat(adatas, join="outer", label=None, index_unique=None)
    validate_unique_obs_names(merged, "merged tissue object")
    merged = merged[:, common_genes].copy()
    merged.var["mt"] = merged.var_names.str.upper().str.startswith("MT-")
    merged.var = merged.var[["mt"]]
    values = merged.X.data if issparse(merged.X) else np.asarray(merged.X)
    if not np.allclose(values, np.round(values)):
        raise ValueError("Merged X is not integer-like; expected raw counts")
    merged.layers["counts"] = merged.X.copy()
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)
    merged.raw = merged.copy()

    output = prepare_output(
        paths.dataset_merge / "tissue_dataset_merged_log_normalized.h5ad",
        args.overwrite,
    )
    merged.write_h5ad(output, compression="gzip")
    print(
        f"Saved {merged.n_obs:,} cells x {merged.n_vars:,} genes to {output}; "
        f"genes required in >= {required}/{len(adatas)} datasets"
    )


if __name__ == "__main__":
    main()
