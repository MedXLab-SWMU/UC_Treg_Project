#!/usr/bin/env python3
"""Stage 5: annotate Immune, Stromal, Epithelial, and Other L1 lineages."""

from __future__ import annotations

import argparse
import re

from uctreg.common import (
    cluster_expression_scores,
    get_cluster_markers,
    run_gpt_celltype,
    safe_genes,
    save_cell_table,
    save_umap,
    set_seed,
    validate_unique_obs_names,
)
from uctreg.configs import L1_CORE, L1_EXT
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root


def normalize_gpt_lineage(label: str) -> str:
    text = re.sub(r"[^a-z0-9 ]", " ", str(label).lower())
    keywords = {
        "Immune": ["immune", "lymphocyte", "t cell", "b cell", "nk cell", "plasma", "macrophage", "monocyte", "dendritic", "myeloid"],
        "Epithelial": ["epithelial", "epithelium", "enterocyte", "goblet", "tuft", "enteroendocrine", "best4"],
        "Stromal": ["stromal", "mesenchymal", "fibroblast", "myofibroblast", "smooth muscle", "endothelial", "pericyte", "perivascular"],
        "Other": ["neural", "neuron", "glia", "schwann", "enteric", "unknown", "uncertain", "unassigned"],
    }
    for lineage, terms in keywords.items():
        if any(term in text for term in terms):
            return lineage
    return "Other"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--run-gpt", action="store_true")
    args = parser.parse_args()
    set_seed(args.seed)

    paths = ProjectPaths(resolve_data_root(args.data_root))
    source = require_input(paths.scvi_integrated / "tissue_dataset_scvi.h5ad")
    output = prepare_output(paths.l1() / "L1_lineage_annotated.h5ad", args.overwrite)
    import pandas as pd
    import scanpy as sc

    adata = sc.read_h5ad(source)
    validate_unique_obs_names(adata, "L1 input")
    cluster_key = "leiden_scVI_0.2"
    if cluster_key not in adata.obs:
        raise KeyError(f"Missing obs['{cluster_key}']")
    adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
    clusters = [str(item) for item in adata.obs[cluster_key].cat.categories]

    panels = {
        lineage: sorted(set(safe_genes(adata, L1_CORE[lineage] + L1_EXT[lineage])))
        for lineage in L1_CORE
    }
    if len(panels["Other"]) < 4:
        del panels["Other"]
    marker_call, marker_gap, _ = cluster_expression_scores(adata, cluster_key, panels)
    adata.obs["L1lin_marker"] = adata.obs[cluster_key].astype(str).map(marker_call).astype("category")
    adata.obs["L1lin_marker_gap"] = adata.obs[cluster_key].astype(str).map(marker_gap).astype(float)

    for lineage, genes in panels.items():
        sc.tl.score_genes(adata, genes, score_name=f"score_L1lin_{lineage}", use_raw=False)
    score_call, score_gap = {}, {}
    for cluster in clusters:
        mask = adata.obs[cluster_key].astype(str) == cluster
        values = pd.Series({lineage: float(adata.obs.loc[mask, f"score_L1lin_{lineage}"].mean()) for lineage in panels}).sort_values(ascending=False)
        score_call[cluster] = str(values.index[0])
        score_gap[cluster] = float(values.iloc[0] - values.iloc[1])
    adata.obs["L1lin_score"] = adata.obs[cluster_key].astype(str).map(score_call).astype("category")
    adata.obs["L1lin_score_gap"] = adata.obs[cluster_key].astype(str).map(score_gap).astype(float)

    gpt_call = {cluster: "Other" for cluster in clusters}
    gpt_raw = {cluster: "not_run" for cluster in clusters}
    if args.run_gpt:
        markers = get_cluster_markers(adata, cluster_key)
        raw = run_gpt_celltype(markers, "colon tissue")
        gpt_raw = {cluster: str(raw.get(cluster, "")) for cluster in clusters}
        gpt_call = {cluster: normalize_gpt_lineage(gpt_raw[cluster]) for cluster in clusters}
    adata.obs["L1lin_gpt"] = adata.obs[cluster_key].astype(str).map(gpt_call).astype("category")

    final_call, confidence = {}, {}
    for cluster in clusters:
        votes = [marker_call[cluster], score_call[cluster]]
        if args.run_gpt:
            votes.append(gpt_call[cluster])
        counts = pd.Series(votes).value_counts()
        final_call[cluster] = str(counts.index[0])
        confidence[cluster] = float(counts.iloc[0] / len(votes))
        if confidence[cluster] < 0.6:
            final_call[cluster] = "Other"
    adata.obs["L1_lineage"] = adata.obs[cluster_key].astype(str).map(final_call).astype("category")
    adata.obs["L1lin_conf"] = adata.obs[cluster_key].astype(str).map(confidence).astype(float)

    cluster_summary = pd.DataFrame({
        "cluster": clusters,
        "n_cells": [int((adata.obs[cluster_key].astype(str) == cluster).sum()) for cluster in clusters],
        "marker": [marker_call[c] for c in clusters],
        "score": [score_call[c] for c in clusters],
        "gpt_raw": [gpt_raw[c] for c in clusters],
        "gpt": [gpt_call[c] for c in clusters],
        "L1_lineage": [final_call[c] for c in clusters],
        "confidence": [confidence[c] for c in clusters],
    })
    paths.l1().mkdir(parents=True, exist_ok=True)
    cluster_summary.to_csv(paths.l1() / "L1_cluster_summary.csv", index=False)
    save_cell_table(
        adata,
        paths.l1() / "L1_lineage_summary.csv",
        ["GSM","group","GSE",cluster_key,"L1lin_marker","L1lin_score","L1lin_gpt","L1_lineage","L1lin_conf"],
    )
    save_umap(adata, paths.l1() / "figs" / "L1_lineage_umap.png", "L1_lineage", title="L1 Lineage Annotation")

    keep = ["GSM","group","GSE",cluster_key,"L1lin_marker","L1lin_score","L1lin_gpt","L1_lineage","L1lin_conf"]
    adata.obs = adata.obs[[column for column in keep if column in adata.obs]].copy()
    adata.obsp.pop("distances", None)
    adata.write_h5ad(output, compression="gzip")
    print(f"Saved L1 annotation to {output}")


if __name__ == "__main__":
    main()
