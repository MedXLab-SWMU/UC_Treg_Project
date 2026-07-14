from __future__ import annotations

import os
import random
from pathlib import Path
from typing import Iterable


def set_seed(seed: int) -> None:
    random.seed(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass


def safe_genes(adata, genes: Iterable[str]) -> list[str]:
    return [gene for gene in genes if gene in adata.var_names]


def require_obs(adata, columns: Iterable[str]) -> None:
    missing = [column for column in columns if column not in adata.obs]
    if missing:
        raise KeyError(f"Missing obs columns: {missing}")


def validate_unique_obs_names(adata, label: str) -> None:
    if not adata.obs_names.is_unique:
        duplicated = adata.obs_names[adata.obs_names.duplicated()].unique()[:5].tolist()
        raise ValueError(f"{label} has duplicate cell barcodes, examples: {duplicated}")


def reset_analysis_state(adata) -> None:
    for key in ("X_pca", "X_pca_harmony", "X_scVI", "X_umap"):
        adata.obsm.pop(key, None)
    adata.uns.pop("neighbors", None)
    for key in ("distances", "connectivities"):
        adata.obsp.pop(key, None)
    for key in ("highly_variable", "highly_variable_rank", "means", "variances", "variances_norm", "highly_variable_nbatches"):
        if key in adata.var:
            del adata.var[key]
    for key in list(adata.uns):
        if "colors" in key or "leiden" in key or key == "rank_genes_groups":
            del adata.uns[key]


def prepare_harmony_subcluster(
    adata,
    *,
    hvg_n: int,
    resolutions: Iterable[float],
    key_suffix: str,
    seed: int,
    n_neighbors: int = 20,
) -> None:
    import scanpy as sc
    import scanpy.external as sce

    reset_analysis_state(adata)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(hvg_n, adata.n_vars),
        flavor="seurat_v3",
        subset=False,
        span=0.5,
        check_values=False,
    )
    sub_hv = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(sub_hv, max_value=10)
    max_components = min(sub_hv.n_vars, sub_hv.n_obs) - 1
    if max_components < 1:
        raise ValueError("At least two cells and two HVGs are required for PCA")
    n_comps = min(50, max_components)
    sc.tl.pca(sub_hv, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    adata.obsm["X_pca"] = sub_hv.obsm["X_pca"].copy()
    sce.pp.harmony_integrate(adata, key="GSE", basis="X_pca")
    sc.pp.neighbors(
        adata,
        use_rep="X_pca_harmony",
        n_neighbors=min(n_neighbors, adata.n_obs - 1),
        random_state=seed,
    )
    sc.tl.umap(adata, random_state=seed)
    for resolution in resolutions:
        key = f"leiden_{resolution}"
        if key_suffix:
            key += f"_{key_suffix}"
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key,
            random_state=seed,
        )


def cluster_expression_scores(adata, group_key: str, panels: dict[str, list[str]]):
    import numpy as np
    import pandas as pd
    from scipy.sparse import issparse

    filtered = {name: safe_genes(adata, genes) for name, genes in panels.items()}
    calls: dict[str, str] = {}
    gaps: dict[str, float] = {}
    categories = adata.obs[group_key].astype("category").cat.categories
    for cluster in categories:
        mask = (adata.obs[group_key] == cluster).to_numpy()
        scores = {}
        for name, genes in filtered.items():
            if not genes:
                scores[name] = float("-inf")
                continue
            matrix = adata[mask, genes].X
            mean = float(matrix.mean())
            positive = float((matrix > 0).mean()) if not issparse(matrix) else float(matrix.getnnz() / (matrix.shape[0] * matrix.shape[1]))
            scores[name] = 0.6 * mean + 0.4 * positive
        ranked = pd.Series(scores).sort_values(ascending=False)
        calls[str(cluster)] = str(ranked.index[0])
        gaps[str(cluster)] = float(ranked.iloc[0] - ranked.iloc[1]) if len(ranked) > 1 else float(ranked.iloc[0])
    return calls, gaps, filtered


def filter_calls_by_gap(calls: dict[str, str], gaps: dict[str, float], quantile: float = 0.2):
    import numpy as np

    values = np.asarray([value for value in gaps.values() if np.isfinite(value)])
    threshold = float(np.quantile(values, quantile)) if values.size else 0.0
    filtered = {cluster: (label if gaps.get(cluster, float("-inf")) >= threshold else "Other") for cluster, label in calls.items()}
    return filtered, threshold


def score_gene_sets(adata, group_key: str, panels: dict[str, list[str]]):
    import pandas as pd
    import scanpy as sc

    for name, genes in panels.items():
        valid = safe_genes(adata, genes)
        if valid:
            sc.tl.score_genes(adata, valid, score_name=f"L2score_{name}", use_raw=False)
    calls: dict[str, str] = {}
    gaps: dict[str, float] = {}
    categories = adata.obs[group_key].astype("category").cat.categories
    for cluster in categories:
        mask = adata.obs[group_key] == cluster
        values = {
            name: float(adata.obs.loc[mask, f"L2score_{name}"].mean())
            if f"L2score_{name}" in adata.obs else float("-inf")
            for name in panels
        }
        ranked = pd.Series(values).sort_values(ascending=False)
        calls[str(cluster)] = str(ranked.index[0])
        gaps[str(cluster)] = float(ranked.iloc[0] - ranked.iloc[1]) if len(ranked) > 1 else float(ranked.iloc[0])
    return calls, gaps


def get_cluster_markers(adata, group_key: str, top_n: int = 10):
    import omicverse as ov
    import scanpy as sc

    sc.tl.rank_genes_groups(adata, groupby=group_key, method="wilcoxon", n_genes=200)
    return ov.single.get_celltype_marker(
        adata,
        clustertype=group_key,
        rank=True,
        key="rank_genes_groups",
        foldchange=1,
        topgenenumber=top_n,
    )


def run_gpt_celltype(markers, tissue_name: str):
    if not os.environ.get("AGI_API_KEY"):
        raise RuntimeError("--run-gpt requires AGI_API_KEY in the environment.")
    import omicverse as ov

    return ov.single.gptcelltype(
        markers,
        tissuename=tissue_name,
        speciename="human",
        model="qwen-plus",
        provider="qwen",
        topgenenumber=10,
    )


def save_cell_table(adata, path: Path, columns: Iterable[str]) -> None:
    columns = [column for column in columns if column in adata.obs]
    table = adata.obs[columns].copy()
    table.index.name = "cell_id"
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=True)


def save_umap(adata, path: Path, color, *, title: str | None = None) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    path.parent.mkdir(parents=True, exist_ok=True)
    sc.pl.umap(adata, color=color, frameon=False, title=title, show=False)
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close("all")
