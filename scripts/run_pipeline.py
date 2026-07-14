#!/usr/bin/env python3
"""Run all or a selected contiguous range of UC Treg pipeline stages."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


STAGES = {
    1: ("01_preprocess_samples.py", ["--dataset", "all"]),
    2: ("02_qc_datasets.py", ["--dataset", "all"]),
    3: ("03_merge_tissue.py", []),
    4: ("04_integrate_scvi.py", []),
    5: ("05_annotate_l1.py", []),
    6: ("06_annotate_l2.py", ["--lineage", "all"]),
    7: ("07_refine_l3_parents.py", ["--lineage", "all"]),
    8: ("08_annotate_l3_immune.py", ["--cell-type", "all"]),
    9: ("09_annotate_l3_stromal.py", ["--cell-type", "all"]),
    10: ("10_finalize_annotations.py", []),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", help="Falls back to UC_TREG_DATA_ROOT.")
    parser.add_argument("--raw-root", help="Falls back to UC_TREG_RAW_ROOT.")
    parser.add_argument("--from-stage", type=int, choices=STAGES, default=1)
    parser.add_argument("--to-stage", type=int, choices=STAGES, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--run-gpt", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.from_stage > args.to_stage:
        parser.error("--from-stage must not exceed --to-stage")
    data_root = args.data_root or os.environ.get("UC_TREG_DATA_ROOT")
    if not data_root:
        parser.error("Pass --data-root or set UC_TREG_DATA_ROOT")

    script_dir = Path(__file__).resolve().parent
    for stage in range(args.from_stage, args.to_stage + 1):
        filename, extras = STAGES[stage]
        command = [sys.executable, str(script_dir / filename), "--data-root", data_root, "--seed", str(args.seed), *extras]
        if stage == 1 and args.raw_root:
            command += ["--raw-root", args.raw_root]
        if args.overwrite:
            command.append("--overwrite")
        if args.run_gpt and stage in {5, 6}:
            command.append("--run-gpt")
        print(f"[stage {stage}] {' '.join(command)}", flush=True)
        if not args.dry_run:
            subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
