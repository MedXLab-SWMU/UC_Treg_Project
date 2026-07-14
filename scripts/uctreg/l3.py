from __future__ import annotations

from pathlib import Path

from .common import prepare_harmony_subcluster, save_cell_table, save_umap, validate_unique_obs_names
from .paths import prepare_output, require_input


def run_l3_branch(
    *,
    source: Path,
    output_dir: Path,
    config: dict,
    overwrite: bool,
    seed: int,
) -> None:
    source = require_input(source)
    stem = config["stem"]
    output_h5ad = prepare_output(output_dir / f"L3_{stem}_annotated.h5ad", overwrite)
    output_csv = prepare_output(output_dir / f"L3_{stem}_summary.csv", overwrite)
    import pandas as pd
    import scanpy as sc

    adata = sc.read_h5ad(source)
    if "L2_celltype_sub" not in adata.obs:
        raise KeyError("Missing obs['L2_celltype_sub']")
    sub = adata[adata.obs["L2_celltype_sub"].astype(str) == config["subset"]].copy()
    if sub.n_obs == 0:
        raise ValueError(f"No cells found for L2 subtype {config['subset']}")
    validate_unique_obs_names(sub, f"L3 {stem}")
    prepare_harmony_subcluster(
        sub,
        hvg_n=config["hvg"],
        resolutions=config["resolutions"],
        key_suffix=config.get("cluster_suffix", stem),
        seed=seed,
    )
    cluster_key = config["cluster_key"]
    if cluster_key not in sub.obs:
        raise KeyError(f"Expected cluster key {cluster_key}; available: {list(sub.obs.columns)}")
    sub.obs[cluster_key] = sub.obs[cluster_key].astype("category")

    sc.tl.rank_genes_groups(sub, groupby=cluster_key, method="wilcoxon", use_raw=True, pts=True)
    top_markers = {}
    for cluster in sub.obs[cluster_key].cat.categories:
        frame = sc.get.rank_genes_groups_df(sub, group=cluster)
        top_markers[str(cluster)] = frame.loc[frame["logfoldchanges"] > 0, "names"].head(10).tolist()
    pd.DataFrame(dict([(key, pd.Series(value)) for key, value in top_markers.items()])).to_csv(
        output_dir / f"L3_{stem}_top10_markers.csv", index=False
    )

    mapping = config.get("mapping")
    if mapping is None:
        sub.obs["L3"] = config["subset"]
    else:
        observed = set(sub.obs[cluster_key].astype(str).unique())
        missing = sorted(observed - set(mapping))
        if missing:
            raise ValueError(
                f"L3 {stem} mapping does not cover clusters {missing}; review markers before continuing"
            )
        sub.obs["L3"] = sub.obs[cluster_key].astype(str).map(mapping)
    sub.obs["L3"] = sub.obs["L3"].astype("category")

    save_cell_table(sub, output_csv, ["GSM", "GSE", "group", "L1_lineage", "L2_celltype_sub", cluster_key, "L3"])
    save_umap(sub, output_dir / "figs" / f"L3_{stem}_umap.png", "L3", title=f"L3 {stem} Annotation")
    sub.write_h5ad(output_h5ad, compression="gzip")
    print(f"Saved L3 {stem}: {sub.n_obs:,} cells -> {output_h5ad}")
