#!/usr/bin/env python3
"""Run the optional fast checks for experiment workers; no benchmark sweep."""
import argparse
from pathlib import Path

from . import ROOT
from .check_optimization_benchmark import check as check_optimizations
from .check_pipeline import check as check_runtime


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["all", "runtime", "optimizations"], default="all")
    parser.add_argument("--runtime-worker", type=Path, default=ROOT / "build/experiments/pruning_experiment")
    parser.add_argument("--optimization-worker", type=Path,
                        default=ROOT / "build/local-release/pruning_optimization_worker")
    args = parser.parse_args()
    if args.suite in ("all", "runtime"):
        check_runtime(args.runtime_worker.resolve())
    if args.suite in ("all", "optimizations"):
        check_optimizations(args.optimization_worker.resolve())


if __name__ == "__main__":
    main()

