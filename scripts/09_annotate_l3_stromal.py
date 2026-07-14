#!/usr/bin/env python3
"""Stage 9: annotate one or all stromal L3 branches."""

from __future__ import annotations

import argparse

from uctreg.common import set_seed
from uctreg.configs import L3_STROMAL
from uctreg.l3 import run_l3_branch
from uctreg.paths import ProjectPaths, add_common_args, resolve_data_root


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_args(parser)
    parser.add_argument("--cell-type", choices=["all", *L3_STROMAL], default="all")
    args = parser.parse_args()
    set_seed(args.seed)
    paths = ProjectPaths(resolve_data_root(args.data_root))
    selected = list(L3_STROMAL) if args.cell_type == "all" else [args.cell_type]
    source = paths.l3("stromal") / "L3_stromal_annotated.h5ad"
    for cell_type in selected:
        run_l3_branch(
            source=source,
            output_dir=paths.l3("stromal"),
            config=L3_STROMAL[cell_type],
            overwrite=args.overwrite,
            seed=args.seed,
        )


if __name__ == "__main__":
    main()
