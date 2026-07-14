#!/usr/bin/env python3
"""Stage 10: merge all L3 branches, validate barcodes, and save the final atlas."""

from __future__ import annotations

import argparse
from pathlib import Path

from uctreg.common import save_umap, set_seed, validate_unique_obs_names
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root


IMMUNE_FILES = [
    "L3_TNK_summary.csv", "L3_B_summary.csv", "L3_Mast_summary.csv",
    "L3_Myeloid_summary.csv", "L3_Plasma_summary.csv",
]
STROMAL_FILES = ["L3_endo_summary.csv", "L3_fibroblast_summary.csv", "L3_SMC_summary.csv"]


def read_annotation_table(path: Path):
    require_input(path, "annotation table")
    import pandas as pd

    frame = pd.read_csv(path)
    if "cell_id" not in frame.columns:
        unnamed = [column for column in frame.columns if column.startswith("Unnamed")]
        if len(unnamed) == 1:
            frame.rename(columns={unnamed[0]: "cell_id"}, inplace=True)
        else:
            raise ValueError(f"{path} has no cell_id column")
    if "L3" not in frame.columns:
        raise ValueError(f"{path} has no L3 column")
    if frame["cell_id"].duplicated().any():
        raise ValueError(f"Duplicate cell_id values within {path}")
    if frame["L3"].isna().any():
        raise ValueError(f"Missing L3 labels in {path}")
    return frame.set_index("cell_id")


def combine_tables(paths: list[Path], label: str):
    frames = [read_annotation_table(path) for path in paths]
    import pandas as pd

    combined = pd.concat(frames, axis=0)
    if combined.index.duplicated().any():
        examples = combined.index[combined.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Overlapping {label} branch barcodes: {examples}")
    return combined


def validate_coverage(full_index, annotation_index) -> None:
    full = set(full_index)
    annotated = set(annotation_index)
    missing = full - annotated
    extra = annotated - full
    if missing or extra:
        raise ValueError(
            f"Annotation coverage mismatch: missing={len(missing)}, extra={len(extra)}; "
            f"missing examples={list(missing)[:5]}, extra examples={list(extra)[:5]}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument(
        "--validate-baseline",
        action="store_true",
        help="Require the cell counts recorded in the original notebooks.",
    )
    args = parser.parse_args()
    set_seed(args.seed)

    paths = ProjectPaths(resolve_data_root(args.data_root))
    merge_dir = paths.merged_annotation()
    outputs = {
        "immune": merge_dir / "immune.annotation.csv",
        "stromal": merge_dir / "stromal.annotation.csv",
        "epithelial": merge_dir / "epithelial.annotation.csv",
        "other": merge_dir / "other.annotation.csv",
        "all": merge_dir / "L3_annotation.csv",
        "h5ad": merge_dir / "L3_annotated.h5ad",
    }
    for output in outputs.values():
        prepare_output(output, args.overwrite)

    # Validate all upstream paths before importing the analysis stack.
    for name in IMMUNE_FILES:
        require_input(paths.l3("immune") / name, "immune annotation table")
    for name in STROMAL_FILES:
        require_input(paths.l3("stromal") / name, "stromal annotation table")
    require_input(paths.l3("epithelial") / "L3_epithelial_summary.csv", "epithelial annotation table")
    require_input(paths.l3("other") / "L3_other_summary.csv", "other annotation table")
    source = require_input(paths.l1() / "L1_lineage_annotated.h5ad")

    import pandas as pd
    import scanpy as sc

    immune = combine_tables([paths.l3("immune") / name for name in IMMUNE_FILES], "immune")
    stromal = combine_tables([paths.l3("stromal") / name for name in STROMAL_FILES], "stromal")
    epithelial = read_annotation_table(paths.l3("epithelial") / "L3_epithelial_summary.csv")
    other = read_annotation_table(paths.l3("other") / "L3_other_summary.csv")
    branches = {"immune": immune, "stromal": stromal, "epithelial": epithelial, "other": other}
    for name, frame in branches.items():
        frame.to_csv(outputs[name], index=True, index_label="cell_id")

    annotation = pd.concat(list(branches.values()), axis=0)
    if annotation.index.duplicated().any():
        examples = annotation.index[annotation.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Overlapping lineage barcodes: {examples}")

    adata = sc.read_h5ad(source)
    validate_unique_obs_names(adata, "final L1 object")
    validate_coverage(adata.obs_names, annotation.index)
    annotation = annotation.loc[adata.obs_names]
    annotation[[column for column in ("GSM", "GSE", "group", "L3") if column in annotation]].to_csv(
        outputs["all"], index=True, index_label="cell_id"
    )

    adata.obs["L3"] = annotation["L3"].astype(str).values
    adata.obs.drop(columns=["L1lin_marker", "L1lin_score", "L1lin_gpt", "L1lin_conf"], errors="ignore", inplace=True)
    adata.uns.pop("rank_genes_groups", None)
    before = adata.n_obs
    counts_before = adata.obs["L3"].value_counts()
    adata = adata[~adata.obs["L3"].isin(["Contaminants", "Low Quality"])].copy()
    adata.obs["L3"] = adata.obs["L3"].astype("category")

    if args.validate_baseline:
        expected = {
            "input_cells": 374507, "contaminants": 40798,
            "low_quality": 21789, "final_cells": 311920,
            "cell_types": 44, "treg": 7340,
        }
        observed = {
            "input_cells": before,
            "contaminants": int(counts_before.get("Contaminants", 0)),
            "low_quality": int(counts_before.get("Low Quality", 0)),
            "final_cells": adata.n_obs,
            "cell_types": adata.obs["L3"].nunique(),
            "treg": int((adata.obs["L3"] == "Treg").sum()),
        }
        if observed != expected:
            raise ValueError(f"Notebook baseline mismatch: expected={expected}, observed={observed}")

    adata.write_h5ad(outputs["h5ad"], compression="gzip")
    save_umap(adata, merge_dir / "figs" / "UC_Colon_CellTypes.png", "L3", title="Cell Types in UC Colon Tissue")
    print(
        f"Final atlas: {before:,} input cells, {adata.n_obs:,} retained, "
        f"{adata.obs['L3'].nunique()} cell types -> {outputs['h5ad']}"
    )


if __name__ == "__main__":
    main()
