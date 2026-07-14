from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path


def resolve_data_root(value: str | os.PathLike[str] | None) -> Path:
    raw = value or os.environ.get("UC_TREG_DATA_ROOT")
    if not raw:
        raise ValueError(
            "Data root is required. Pass --data-root or set UC_TREG_DATA_ROOT."
        )
    return Path(raw).expanduser().resolve()


def resolve_raw_root(
    value: str | os.PathLike[str] | None, data_root: Path
) -> Path:
    raw = value or os.environ.get("UC_TREG_RAW_ROOT")
    return Path(raw).expanduser().resolve() if raw else data_root / "01_raw"


@dataclass(frozen=True)
class ProjectPaths:
    root: Path

    @property
    def sample_merge(self) -> Path:
        return self.root / "02_sample_merge"

    @property
    def qc_filtered(self) -> Path:
        return self.root / "03_qc_filtered"

    @property
    def dataset_merge(self) -> Path:
        return self.root / "04_dataset_merge"

    @property
    def scvi_integrated(self) -> Path:
        return self.root / "05_scvi_integrated"

    @property
    def annotation(self) -> Path:
        return self.root / "06_annotation"

    def l1(self) -> Path:
        return self.annotation / "L1"

    def l2(self) -> Path:
        return self.annotation / "L2"

    def l3(self, lineage: str) -> Path:
        return self.annotation / "L3" / lineage

    def merged_annotation(self) -> Path:
        return self.annotation / "Merge"


def require_input(path: Path, label: str = "input") -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


def prepare_output(path: Path, overwrite: bool) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}. Pass --overwrite to replace it."
        )
    return path


def add_common_args(parser: argparse.ArgumentParser, *, raw_root: bool = False) -> None:
    parser.add_argument(
        "--data-root",
        help="Pipeline data root; falls back to UC_TREG_DATA_ROOT.",
    )
    if raw_root:
        parser.add_argument(
            "--raw-root",
            help="Raw-data root; falls back to UC_TREG_RAW_ROOT or <data-root>/01_raw.",
        )
    parser.add_argument("--overwrite", action="store_true", help="Replace outputs.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed.")
