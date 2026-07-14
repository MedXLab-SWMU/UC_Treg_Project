#!/usr/bin/env python3
"""Stage 4: train scVI on batch-aware HVGs and build the integrated graph."""

from __future__ import annotations

import argparse

from uctreg.common import require_obs, set_seed, validate_unique_obs_names
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--hvg", type=int, default=3000)
    parser.add_argument("--latent", type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.2)
    args = parser.parse_args()
    set_seed(args.seed)

    paths = ProjectPaths(resolve_data_root(args.data_root))
    source = require_input(paths.dataset_merge / "tissue_dataset_merged_log_normalized.h5ad")
    output = prepare_output(paths.scvi_integrated / "tissue_dataset_scvi.h5ad", args.overwrite)
    import numpy as np
    import omicverse as ov
    import scanpy as sc
    import scvi

    scvi.settings.seed = args.seed
    adata = sc.read_h5ad(source)
    validate_unique_obs_names(adata, "merged tissue object")
    require_obs(adata, ["GSE", "GSM", "group"])
    if "counts" not in adata.layers:
        raise KeyError("Missing adata.layers['counts']")
    adata.obs = adata.obs[["GSM", "group", "GSE"]].copy()
    for column in ("GSE", "GSM"):
        adata.obs[column] = adata.obs[column].astype("category")

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(args.hvg, adata.n_vars),
        flavor="seurat_v3",
        batch_key="GSE",
        layer="counts",
        inplace=True,
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    ov.single.batch_correction(
        adata_hvg,
        batch_key="GSM",
        methods="scVI",
        n_layers=2,
        n_latent=args.latent,
        gene_likelihood="nb",
    )
    if "X_scVI" not in adata_hvg.obsm:
        raise RuntimeError("scVI did not produce adata_hvg.obsm['X_scVI']")
    if not np.array_equal(adata.obs_names, adata_hvg.obs_names):
        raise ValueError("Cell order changed while creating the HVG object")
    adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"].copy()

    for column in ("highly_variable", "highly_variable_rank", "means", "variances", "variances_norm", "highly_variable_nbatches"):
        if column in adata.var:
            del adata.var[column]
    adata.uns.pop("hvg", None)
    sc.pp.neighbors(adata, use_rep="X_scVI", random_state=args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    cluster_key = f"leiden_scVI_{args.resolution}"
    sc.tl.leiden(adata, resolution=args.resolution, key_added=cluster_key, random_state=args.seed)
    adata.write_h5ad(output, compression="gzip")
    print(f"Saved integrated object {adata.shape} to {output}")


if __name__ == "__main__":
    main()
