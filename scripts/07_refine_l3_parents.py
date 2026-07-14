#!/usr/bin/env python3
"""Stage 7: resolve L2 Other clusters before lineage-specific L3 analysis."""

from __future__ import annotations

import argparse

from uctreg.common import save_cell_table, save_umap, set_seed, validate_unique_obs_names
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root


LINEAGES = ("immune", "stromal", "epithelial", "other")


def refine(lineage: str, paths: ProjectPaths, overwrite: bool) -> None:
    source = require_input(paths.l2() / f"L2_{lineage}_annotated.h5ad")
    output_dir = paths.l3(lineage)
    output_h5ad = prepare_output(output_dir / f"L3_{lineage}_annotated.h5ad", overwrite)
    output_csv = prepare_output(output_dir / f"L3_{lineage}_summary.csv", overwrite)
    import scanpy as sc

    adata = sc.read_h5ad(source)
    validate_unique_obs_names(adata, f"L3 parent {lineage}")
    labels = adata.obs["L2_celltype_sub"].astype(str)
    if lineage == "immune":
        if "leiden_0.4" not in adata.obs:
            raise KeyError("Immune refinement requires obs['leiden_0.4']")
        clusters = adata.obs["leiden_0.4"].astype(str)
        labels.loc[clusters.isin(["5", "7", "10"])] = "T_NK"
        labels.loc[clusters == "8"] = "Mast"
    elif lineage == "stromal":
        labels = labels.replace({"Other": "Fibroblast"})
    elif lineage == "epithelial":
        labels = labels.replace({"Other": "Early_Absorptive"})
    adata.obs["L2_celltype_sub"] = labels.astype("category")
    adata.obs["L3"] = adata.obs["L2_celltype_sub"].astype(str).astype("category")
    save_cell_table(adata, output_csv, ["GSM", "GSE", "group", "L1_lineage", "L2_celltype_sub", "L3"])
    save_umap(adata, output_dir / "figs" / f"L3_{lineage}_parent_umap.png", "L2_celltype_sub")
    adata.write_h5ad(output_h5ad, compression="gzip")
    print(f"Saved refined {lineage} parent object -> {output_h5ad}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--lineage", choices=["all", *LINEAGES], default="all")
    args = parser.parse_args()
    set_seed(args.seed)
    paths = ProjectPaths(resolve_data_root(args.data_root))
    selected = LINEAGES if args.lineage == "all" else (args.lineage,)
    for lineage in selected:
        refine(lineage, paths, args.overwrite)


if __name__ == "__main__":
    main()
