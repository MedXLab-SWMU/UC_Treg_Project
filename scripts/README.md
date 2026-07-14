# UC Treg scRNA-seq scripts

These scripts are the compact, executable form of the notebooks. The notebooks remain the analysis record; the scripts are the maintained pipeline.

## Data layout

Pass the external data root with `--data-root` or set `UC_TREG_DATA_ROOT`. Raw inputs default to `<data-root>/01_raw`; use `--raw-root` or `UC_TREG_RAW_ROOT` when raw data live elsewhere.

The scripts preserve the existing stage directories:

```text
01_raw/                 optional default raw-data location
02_sample_merge/        merged samples after Scrublet
03_qc_filtered/         strict dataset-level QC
04_dataset_merge/       merged colon counts/log-normalized data
05_scvi_integrated/     scVI latent space and global clustering
06_annotation/          L1, L2, L3, and final merged annotations
```

## Running

Run the full pipeline:

```bash
python scripts/run_pipeline.py --data-root /path/to/scRNA-seq_data --raw-root /path/to/raw
```

Run a range or preview commands:

```bash
python scripts/run_pipeline.py --data-root /path/to/data --from-stage 3 --to-stage 6
python scripts/run_pipeline.py --data-root /path/to/data --dry-run
```

Every stage can run independently. Examples:

```bash
python scripts/01_preprocess_samples.py --data-root /path/to/data --raw-root /path/to/raw --dataset 1_GSE182270
python scripts/06_annotate_l2.py --data-root /path/to/data --lineage immune
python scripts/08_annotate_l3_immune.py --data-root /path/to/data --cell-type tnk
```

Existing outputs are protected. Pass `--overwrite` only when replacement is intended.

GPT annotation is disabled by default because it is external and billable. To enable it:

```bash
export AGI_API_KEY='...'
python scripts/05_annotate_l1.py --data-root /path/to/data --run-gpt
```

No API key is stored in these scripts.

When GPT is off, L2 uses the cluster decisions already reviewed and saved in
the notebooks, while still recalculating marker and gene-set evidence. If the
new clustering no longer matches those reviewed cluster IDs, the script stops
for manual review instead of applying stale labels.

Use the same Python environment as the notebooks. Required analysis packages are
`scanpy`, `anndata`, `omicverse`, `scvi-tools`, `harmonypy`, `leidenalg`,
`pandas`, `numpy`, `scipy`, and `matplotlib`.

## Important behavior

- Tissue merging excludes PBMC and uses `ceil(n_datasets * 0.8)`. For 11 tissue datasets, a gene must occur in at least 9.
- Branch annotation tables are joined by `cell_id`; duplicates, overlaps, missing cells, and extra cells stop the pipeline.
- L3 mappings are intentionally fixed to those reviewed in the notebooks. A new/unmapped Leiden cluster stops the script for manual marker review.
- Stage 10 can compare the result with the original notebook counts using `--validate-baseline`.
