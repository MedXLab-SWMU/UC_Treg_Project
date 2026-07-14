# UC Treg Project

Integrated single-cell RNA-seq analysis of ulcerative colitis colon datasets.
The original analysis record is retained under `notebooks/`; the compact,
executable pipeline and usage instructions are under [`scripts/`](scripts/README.md).

The scripted workflow covers dataset loading, doublet/QC filtering, tissue
merging, scVI integration, L1-L3 annotation, barcode-validated annotation
merging, and final low-quality/contaminant removal.
