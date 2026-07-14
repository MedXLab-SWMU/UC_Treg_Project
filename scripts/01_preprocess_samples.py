#!/usr/bin/env python3
"""Stage 1: read each public dataset, merge samples, and remove doublets."""

from __future__ import annotations

import argparse
from pathlib import Path

from uctreg.common import set_seed, validate_unique_obs_names
from uctreg.configs import (
    DATASET_ORDER,
    DATASET_SPECS,
    GSE116222_GSM_MAP,
    GSE116222_ID_MAP,
    SCP259_SAMPLE_TO_ID,
)
from uctreg.paths import ProjectPaths, add_common_args, prepare_output, require_input, resolve_data_root, resolve_raw_root


def _read_10x_directory(path: Path):
    import scanpy as sc

    adatas = []
    for sample_dir in sorted(item for item in path.iterdir() if item.is_dir()):
        adata = sc.read_10x_mtx(sample_dir, var_names="gene_symbols", cache=False)
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata.obs["GSM_id"] = sample_dir.name
        adata.obs_names = [f"{sample_dir.name}_{barcode}" for barcode in adata.obs_names]
        adatas.append(adata)
    if not adatas:
        raise FileNotFoundError(f"No 10x sample directories found under {path}")
    adata = sc.concat(adatas, merge="same")
    adata.obs["GSM"] = adata.obs["GSM_id"].str.split("_").str[0]
    adata.obs["id"] = adata.obs["GSM_id"].str.split("_").str[1]
    return adata


def _matrix_files(path: Path):
    return sorted(
        item for item in path.iterdir()
        if item.is_file() and item.suffix.lower() in {".txt", ".tsv", ".csv"}
    )


def _read_matrix_directory(path: Path, *, transpose: bool, prefix_sample: bool = True):
    import anndata as ad
    import pandas as pd
    import scanpy as sc

    adatas = []
    for matrix_path in _matrix_files(path):
        sep = "," if matrix_path.suffix.lower() == ".csv" else "\t"
        frame = pd.read_csv(matrix_path, sep=sep, index_col=0)
        matrix = frame.T if transpose else frame
        sample = matrix_path.stem
        item = ad.AnnData(matrix)
        sc.pp.filter_cells(item, min_genes=200)
        sc.pp.filter_genes(item, min_cells=3)
        if prefix_sample:
            item.obs["GSM_id"] = sample
            item.obs_names = [f"{sample}_{barcode}" for barcode in item.obs_names]
        adatas.append(item)
    if not adatas:
        raise FileNotFoundError(f"No expression matrices found under {path}")
    merged = sc.concat(adatas, merge="same")
    if prefix_sample:
        merged.obs["GSM"] = merged.obs["GSM_id"].str.split("_").str[0]
        merged.obs["id"] = merged.obs["GSM_id"].str.split("_").str[1]
    return merged


def _read_gse116222(path: Path):
    adata = _read_matrix_directory(path, transpose=True, prefix_sample=False)
    adata.obs["source_group"] = adata.obs_names.to_series().str.extract(r"([A-Z]\d+)", expand=False).str.strip().values
    adata.obs["id"] = adata.obs["source_group"].map(GSE116222_ID_MAP)
    adata.obs["GSM"] = adata.obs["id"].map(GSE116222_GSM_MAP)
    return adata


def _read_gse235663(path: Path):
    adata = _read_matrix_directory(path, transpose=True, prefix_sample=False)
    batch_code = adata.obs_names.to_series().str.split(".").str[-1].astype(int).values
    adata.obs["batch_code"] = batch_code
    adata = adata[~adata.obs["batch_code"].isin([1, 4])].copy()
    id_map = {2: "UC75", 3: "HC95", 5: "HC96", 6: "UC77"}
    gsm_map = {2: "GSM7507205", 3: "GSM7507206", 5: "GSM7507217", 6: "GSM7507220"}
    adata.obs["id"] = adata.obs["batch_code"].map(id_map)
    adata.obs["GSM"] = adata.obs["batch_code"].map(gsm_map)
    del adata.obs["batch_code"]
    return adata


def _read_gse114374(path: Path):
    import scanpy as sc

    parts = []
    mappings = {
        "HC": ({"S66": "GSM3140593", "S90": "GSM3140594"}, {"S66": "HC87", "S90": "HC88"}),
        "UC": ({"S78": "GSM3140595", "S54": "GSM3140596"}, {"S78": "UC66", "S54": "UC67"}),
    }
    for group, (gsm_map, id_map) in mappings.items():
        part = _read_matrix_directory(require_input(path / group), transpose=True, prefix_sample=False)
        part.obs["sample"] = part.obs_names.str.split("-").str[-1]
        part.obs["GSM"] = part.obs["sample"].map(gsm_map)
        part.obs["id"] = part.obs["sample"].map(id_map)
        parts.append(part)
    return sc.concat(parts, merge="same")


def _read_scp259(path: Path):
    import scanpy as sc

    adata = sc.read_h5ad(path)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    if "Sample" not in adata.obs:
        raise KeyError("SCP259 input must contain obs['Sample']")
    adata.obs["id"] = adata.obs["Sample"].map(SCP259_SAMPLE_TO_ID)
    if "GSM" not in adata.obs:
        adata.obs["GSM"] = adata.obs["Sample"].astype(str)
    return adata


def load_dataset(dataset: str, raw_path: Path):
    loader = DATASET_SPECS[dataset]["loader"]
    if loader == "10x":
        return _read_10x_directory(raw_path)
    if loader == "matrix_t":
        return _read_matrix_directory(raw_path, transpose=True)
    if loader == "matrix":
        return _read_matrix_directory(raw_path, transpose=False)
    if loader == "gse116222":
        return _read_gse116222(raw_path)
    if loader == "gse235663":
        return _read_gse235663(raw_path)
    if loader == "gse114374":
        return _read_gse114374(raw_path)
    if loader == "scp259":
        return _read_scp259(raw_path)
    raise ValueError(f"Unsupported loader: {loader}")


def process_dataset(dataset: str, raw_root: Path, paths: ProjectPaths, overwrite: bool) -> None:
    raw_path = require_input(raw_root / DATASET_SPECS[dataset]["raw"], "raw dataset")
    output = prepare_output(paths.sample_merge / f"{dataset}.h5ad", overwrite)
    import omicverse as ov

    adata = load_dataset(dataset, raw_path)
    if adata.obs[["GSM", "id"]].isna().any().any():
        raise ValueError(f"{dataset}: incomplete GSM/id metadata mapping")
    validate_unique_obs_names(adata, dataset)
    adata = ov.pp.qc(
        adata,
        tresh={"mito_perc": 0.2, "nUMIs": 500, "detected_genes": 250},
        doublets_method="scrublet",
        batch_key=None,
    )
    adata.write_h5ad(output, compression="gzip")
    print(f"[{dataset}] saved {adata.n_obs:,} cells to {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser, raw_root=True)
    parser.add_argument("--dataset", choices=["all", *DATASET_ORDER], default="all")
    args = parser.parse_args()
    data_root = resolve_data_root(args.data_root)
    raw_root = resolve_raw_root(args.raw_root, data_root)
    paths = ProjectPaths(data_root)
    set_seed(args.seed)
    selected = DATASET_ORDER if args.dataset == "all" else [args.dataset]
    for dataset in selected:
        process_dataset(dataset, raw_root, paths, args.overwrite)


if __name__ == "__main__":
    main()
