#!/usr/bin/env python3
"""Run the paper's size/epsilon experiment: python3 benchmark_runtime.py."""
import argparse

from support import RUNTIME_SIZES
from support.paper_benchmarks import main

# Experiment settings. No command-line arguments are needed.
POINT_COUNTS = RUNTIME_SIZES
EPSILON_PERCENTAGES = [1, 2, 4, 8]
REPEATS = 3
TIMEOUT_SECONDS = 600
MEMORY_MB = 12288

if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__, epilog="Settings are defined at the top of this file.").parse_args()
    main("runtime", points=POINT_COUNTS, epsilon_percents=EPSILON_PERCENTAGES,
         repeats=REPEATS, timeout=TIMEOUT_SECONDS, memory_mb=MEMORY_MB)
