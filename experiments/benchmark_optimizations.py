#!/usr/bin/env python3
"""Run the paper's one-at-a-time experiment: python3 benchmark_optimizations.py."""
import argparse
import subprocess
import sys

from support import HERE
from support.paper_benchmarks import main

# Experiment settings. Each switch is compared with the all-off baseline.
POINT_COUNTS = [80, 160, 320, 640]
EPSILON_PERCENTAGES = [1, 2, 5, 10]
REPEATS = 3
TIMEOUT_SECONDS = 600
MEMORY_MB = 12288
DATA_DIRECTORY = HERE / "optimization_data"

if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__, epilog="Settings are defined at the top of this file.").parse_args()
    if not (DATA_DIRECTORY / "manifest.csv").is_file():
        subprocess.run([sys.executable, str(HERE / "generate_data.py"),
                        "--sizes", *map(str, POINT_COUNTS), "--output-dir", str(DATA_DIRECTORY)], check=True)
    main("optimizations", points=POINT_COUNTS, epsilon_percents=EPSILON_PERCENTAGES,
         repeats=REPEATS, timeout=TIMEOUT_SECONDS, memory_mb=MEMORY_MB,
         manifest=DATA_DIRECTORY / "manifest.csv")
