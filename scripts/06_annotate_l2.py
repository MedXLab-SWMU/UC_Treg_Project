#!/usr/bin/env python3
"""Stage 6: parameterized L2 annotation for one or all L1 lineages."""

from __future__ import annotations

import argparse
import re

from uctreg.common import (
    cluster_expression_scores,
    filter_calls_by_gap,
    get_cluster_markers,
    prepare_harmony_subcluster,
    run_gpt_celltype,
    save_cell_table,
    save_umap,
    score_gene_sets,
    set_seed,
    validate_unique_obs_names,
)
from uctreg.configs import L2_CONFIGS, L2_REVIEWED_CLUSTER_MAPS
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root


def normalize_gpt_label(lineage: str, label: str) -> str:
    text = re.sub(r"[^a-z0-9 ]", " ", str(label).lower())
    terms = {
        "immune": {
            "Plasma": ["plasma"], "Myeloid": ["macroph", "monocyte", "myeloid", "dendritic"],
            "B": ["b cell", "b-cell"], "T_NK": ["t cell", "t-cell", "nk", "lymphocyte"],
        },
        "stromal": {
            "Fibroblast": ["fibro", "mesench"], "Endothelial": ["endothelial"],
            "Pericyte_SMC": ["pericyte", "smooth muscle", "smc", "myofibro"],
        },
        "epithelial": {
            "Absorptive": ["absorptive", "enterocyte", "colonocyte", "best4"],
            "Goblet": ["goblet"], "Enteroendocrine": ["enteroendocrine", "eec"],
            "Tuft": ["tuft"], "Stem_TA_prolif": ["stem", "transit amplifying", "prolifer"],
        },
        "other": {
            "Active_Glia": ["active glia", "activated glia", "stress"],
            "Glia-like_Stromal": ["stromal", "mesenchymal"],
            "Glia": ["glia", "glial", "schwann", "enteric"],
        },
    }
    for category, keywords in terms[lineage].items():
        if any(keyword in text for keyword in keywords):
            return category
    return "Other"


def vote_labels(clusters, marker_call, score_call, gpt_call, run_gpt):
    import pandas as pd

    final, confidence = {}, {}
    for cluster in clusters:
        votes = [marker_call[cluster], score_call[cluster]]
        if run_gpt:
            votes.append(gpt_call[cluster])
        votes = [vote for vote in votes if vote != "Other"]
        if len(votes) < 2:
            final[cluster], confidence[cluster] = "Other", 0.0
            continue
        counts = pd.Series(votes).value_counts()
        confidence[cluster] = float(counts.iloc[0] / len(votes))
        final[cluster] = str(counts.index[0]) if confidence[cluster] >= 2 / 3 else "Other"
    return final, confidence


def annotate_lineage(lineage: str, paths: ProjectPaths, args) -> None:
    config = L2_CONFIGS[lineage]
    source = require_input(paths.l1() / "L1_lineage_annotated.h5ad")
    output = prepare_output(paths.l2() / f"L2_{lineage}_annotated.h5ad", args.overwrite)
    import numpy as np
    import pandas as pd
    import scanpy as sc

    adata = sc.read_h5ad(source)
    sub = adata[adata.obs["L1_lineage"].astype(str) == config.l1_label].copy()
    if sub.n_obs == 0:
        raise ValueError(f"No cells found for L1 lineage {config.l1_label}")
    validate_unique_obs_names(sub, f"L2 {lineage}")
    prepare_harmony_subcluster(
        sub,
        hvg_n=2000,
        resolutions=config.resolutions,
        key_suffix="",
        seed=args.seed,
    )
    main_key = f"leiden_{config.main_resolution}"
    sub.obs[main_key] = sub.obs[main_key].astype("category")
    clusters = [str(item) for item in sub.obs[main_key].cat.categories]

    raw_marker, marker_gap, panels = cluster_expression_scores(sub, main_key, config.panels)
    marker_call, marker_threshold = filter_calls_by_gap(
        raw_marker, marker_gap, config.marker_gap_quantile
    )
    raw_score, score_gap = score_gene_sets(sub, main_key, panels)
    score_call, score_threshold = filter_calls_by_gap(
        raw_score, score_gap, config.score_gap_quantile
    )
    sub.obs["L2_marker"] = sub.obs[main_key].astype(str).map(marker_call).astype("category")
    sub.obs["L2_marker_gap"] = sub.obs[main_key].astype(str).map(marker_gap).astype(float)
    sub.obs["L2_score"] = sub.obs[main_key].astype(str).map(score_call).astype("category")
    sub.obs["L2_score_gap"] = sub.obs[main_key].astype(str).map(score_gap).astype(float)

    gpt_raw = {cluster: "not_run" for cluster in clusters}
    gpt_call = {cluster: "Other" for cluster in clusters}
    if args.run_gpt:
        markers = get_cluster_markers(sub, main_key)
        raw = run_gpt_celltype(markers, config.tissue_name)
        gpt_raw = {cluster: str(raw.get(cluster, "")) for cluster in clusters}
        gpt_call = {cluster: normalize_gpt_label(lineage, gpt_raw[cluster]) for cluster in clusters}
        final, method_conf = vote_labels(clusters, marker_call, score_call, gpt_call, True)
        vote_source = "marker+score+gpt"
    else:
        reviewed = L2_REVIEWED_CLUSTER_MAPS[lineage]
        if set(clusters) != set(reviewed):
            raise ValueError(
                f"L2 {lineage} clusters differ from the reviewed notebook mapping: "
                f"observed={clusters}, reviewed={sorted(reviewed)}. Re-run with --run-gpt "
                "or review and update the mapping."
            )
        final = {cluster: reviewed[cluster] for cluster in clusters}
        method_conf = {
            cluster: float(
                (marker_call[cluster] == final[cluster])
                + (score_call[cluster] == final[cluster])
            ) / 2
            for cluster in clusters
        }
        vote_source = "reviewed_notebook+marker+score"
    sub.obs["L2_gpt"] = sub.obs[main_key].astype(str).map(gpt_call).astype("category")
    sub.obs["L2_celltype_sub"] = pd.Categorical(
        sub.obs[main_key].astype(str).map(final), categories=config.categories
    )
    sub.obs["L2_conf_methods"] = sub.obs[main_key].astype(str).map(method_conf).astype(float)
    sub.obs["L2_votes"] = vote_source
    sub.obs["L2_celltype"] = config.prefix + sub.obs["L2_celltype_sub"].astype(str)

    # Re-score the alternate resolutions with the same two non-external methods.
    resolution_labels = {}
    for resolution in config.resolutions:
        key = f"leiden_{resolution}"
        if key == main_key:
            resolution_labels[resolution] = sub.obs["L2_celltype_sub"].astype(str)
            continue
        alt_marker, _, _ = cluster_expression_scores(sub, key, panels)
        alt_score, _ = score_gene_sets(sub, key, panels)
        alt_map = {
            cluster: alt_marker[cluster] if alt_marker[cluster] == alt_score[cluster] else "Other"
            for cluster in alt_marker
        }
        resolution_labels[resolution] = sub.obs[key].astype(str).map(alt_map)
    stability = {}
    for cluster in clusters:
        mask = sub.obs[main_key].astype(str) == cluster
        target = final[cluster]
        stability[cluster] = float(np.mean([
            (resolution_labels[resolution][mask] == target).mean()
            for resolution in config.resolutions
        ]))
    sub.obs["L2_stability"] = sub.obs[main_key].astype(str).map(stability).astype(float)
    sub.obs["L2_conf"] = 0.8 * sub.obs["L2_conf_methods"] + 0.2 * sub.obs["L2_stability"]

    summary = pd.DataFrame({
        "cluster": clusters,
        "n_cells": [int((sub.obs[main_key].astype(str) == cluster).sum()) for cluster in clusters],
        "marker": [marker_call[c] for c in clusters], "marker_gap": [marker_gap[c] for c in clusters],
        "score": [score_call[c] for c in clusters], "score_gap": [score_gap[c] for c in clusters],
        "gpt_raw": [gpt_raw[c] for c in clusters], "gpt": [gpt_call[c] for c in clusters],
        "L2": [final[c] for c in clusters], "method_conf": [method_conf[c] for c in clusters],
        "stability": [stability[c] for c in clusters],
    })
    paths.l2().mkdir(parents=True, exist_ok=True)
    summary.to_csv(paths.l2() / f"L2_{lineage}_cluster_summary.csv", index=False)
    save_cell_table(
        sub, paths.l2() / f"L2_{lineage}_summary.csv",
        ["GSM","GSE","group","leiden_scVI_0.2","L1_lineage",main_key,"L2_marker","L2_score","L2_gpt","L2_celltype_sub","L2_celltype","L2_conf_methods","L2_stability","L2_conf","L2_votes"],
    )
    save_umap(sub, paths.l2() / "figs" / f"L2_{lineage}_umap.png", "L2_celltype", title=f"L2 {lineage.title()} Annotation")
    keep = ["GSM","GSE","group","leiden_scVI_0.2","L1_lineage",main_key,"L2_marker","L2_score","L2_gpt","L2_celltype_sub","L2_celltype","L2_conf_methods","L2_stability","L2_conf","L2_votes"]
    sub.obs = sub.obs[[column for column in keep if column in sub.obs]].copy()
    sub.write_h5ad(output, compression="gzip")
    print(
        f"Saved L2 {lineage} to {output}; marker gap threshold={marker_threshold:.4g}, "
        f"score gap threshold={score_threshold:.4g}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--lineage", choices=["all", *L2_CONFIGS], default="all")
    parser.add_argument("--run-gpt", action="store_true")
    args = parser.parse_args()
    set_seed(args.seed)
    paths = ProjectPaths(resolve_data_root(args.data_root))
    selected = list(L2_CONFIGS) if args.lineage == "all" else [args.lineage]
    for lineage in selected:
        annotate_lineage(lineage, paths, args)


if __name__ == "__main__":
    main()
