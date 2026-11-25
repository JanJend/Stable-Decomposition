#!/usr/bin/env python3
"""Run the paper's one-at-a-time experiment: python3 benchmark_optimizations.py."""
import argparse

from support.paper_benchmarks import main
from support.generator_options import add_generator_arguments, generator_sweep

# Each optional improvement is compared with a baseline using A-1/2 and optional
# switches off (including the single-solve quotient and updated I target); Full End omits the
# baseline's quotient to measure its saving.
# Dense pivots, fused I, deferred K minimization and throttled progress always run.
POINT_COUNTS = [160, 320, 640]
FAMILIES = ["torus_gauss", "immuno", "gauss_clust"]
EPSILON_PERCENTAGES = [0.5,1,2]
REPEATS = 3
TIMEOUT_SECONDS = 1800
MEMORY_MB = 12288
ORBIT_ROUNDS = [12,24]
ORBIT_SEED_MODES = ["hybrid"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, epilog="Dataset settings are defined at the top of this file.")
    add_generator_arguments(parser, rounds=ORBIT_ROUNDS, modes=ORBIT_SEED_MODES)
    args = parser.parse_args()
    main("optimizations", points=POINT_COUNTS, epsilon_percents=EPSILON_PERCENTAGES,
         families=FAMILIES, repeats=REPEATS, timeout=TIMEOUT_SECONDS, memory_mb=MEMORY_MB,
         generator_configs=generator_sweep(args))
