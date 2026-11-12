# Research experiments

Run either experiment from the repository root, with **no arguments**:

```sh
python3 experiments/benchmark_runtime.py
python3 experiments/benchmark_optimizations.py
```

The runtime experiment uses all ten optimizations. The optimization experiment
compares the all-off baseline with **each of the ten switches enabled alone**.
Both use the existing Release executable, `build/pruning`, with compact columns.
Build it with CMake first; the benchmarks never configure or compile code.

The runtime experiment sweeps target sizes **20, 32, 64, 128, 256, 512, 1024,
2048**, using `experiments/manifest.csv`. The optimization comparison uses only
**80, 160, 320, 640**, because the baseline is expensive on larger inputs.
It prepares those synthetic and observed real inputs in
`experiments/optimization_data/` on the first run and reuses them afterwards.
Both experiments sweep epsilon values **1%, 2%, 5%, 10%** of each presentation's
largest degree extent. Actual sizes, seeds, homology degrees, scalar epsilons
and percentages are retained in the tables. Real clouds use the nearest
available complete observed sample rather than an exact-size subsample.

Results go to `experiments/paper_results/runtime/` and
`experiments/paper_results/optimizations/`, respectively. Each contains
**`table.csv` and `table.tex`**, updated after every measurement. Rerunning the
same script continues an interrupted experiment or regenerates completed tables.
If the code, inputs or settings changed, a timestamped result directory is used
and earlier measurements are preserved. Rerunning also finds and continues that
matching timestamped run. Ctrl-C stops the current measurement and keeps
completed results.

The few adjustable settings—point counts, epsilon percentages, repetitions,
timeout and memory limit—are at the top of each script. No benchmark options
are required. Run the scripts one at a time so their timings do not compete.

## Build and prerequisites

Use Python 3.9 or newer, CMake and the normal C++ project dependencies. Only data
generation needs NumPy (`python3 -m pip install numpy`), including the first
optimization run if its inputs are not yet prepared. Timing and reporting use
the standard library. Build the executable separately before benchmarking:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGRLINA_USE_COMPACT_COLUMNS=ON
cmake --build build --target pruning --config Release -j 2
```

The sibling Persistence-Algebra and AIDA checkouts supply the normal project
dependencies; AIDA decomposition is not run by these benchmarks. On Apple
Silicon, configure with `-DCMAKE_PREFIX_PATH=/opt/homebrew` first if needed.
No second build or vector-storage executable is needed.

Prepare the data once if `experiments/manifest.csv` is absent. For module
generation, build function-Delaunay in its sibling checkout:

```sh
cmake -S ../function_delaunay -B ../function_delaunay/build -DCMAKE_BUILD_TYPE=Release
cmake --build ../function_delaunay/build --target main -j 2
```

Observed immune-cell inputs come from the sibling
`SpatialPatterningOfImmuneCells/1.5mmRegions/` checkout. Its source is
[SpatialPatterningOfImmuneCells](https://github.com/MultiparameterTDAHistology/SpatialPatterningOfImmuneCells).
Pass `--immune-repo PATH` to `generate_data.py` for another location.

## 1. Prepare the datasets

```sh
python3 experiments/generate_data.py
```

Default target sizes are **20, 32, 64, 128, 256, 512, 1024, 2048**.
Every synthetic family/size/seed is a separate deterministic sample.
Real clouds are selected once, independently of the synthetic seeds.

| Family | Sampling | Default homology |
| --- | --- | --- |
| `torus_gaussian` | Area-uniform torus in R3, major radius 3, minor radius 1, Gaussian noise sigma 0.1 | H1 |
| `torus_uniform` | Same noisy torus, replacing 10% of points with ambient-box noise in [-4,4]^2 x [-1,1] | H1 |
| `gaussian_clusters` | Three planar Gaussian clusters; weights 0.4/0.35/0.25, standard deviations 0.3/0.4/0.25 | H0 |
| `sphere4d` | Uniform unit 3-sphere in R4 with Gaussian noise sigma 0.05 | H1 |
| `immuno` | Distinct complete observed immune-cell clouds nearest the target sizes | H1 |

The real data are never subsampled or fabricated. Every CSV under the observed
point-cloud directory is considered, with headers excluded and coincident
points removed. Selection minimizes relative size difference. Actual counts
and observed filenames are retained in the manifest; the largest selected
cloud has 2055 points for target 2048.

Spatial coordinates are translated and isotropically scaled so their largest
coordinate range is one. The final point-set column is negative Gaussian KDE,
linearly normalized to [-1,0], including self-contributions. Bandwidth is mean
spatial standard deviation times n^(-1/(d+4)); exact quadratic KDE uses blocks.
Uniform-noise points participate in the KDE. The R4 input still gives a
two-parameter distance/codensity module, not a four-parameter filtration.

Outputs are `point-sets/<family>/`, `test_presentations/<family>/`, and
`manifest.csv` in the chosen output directory. Function-Delaunay computes
minimal presentations with `--minpres H --no-delaunay-compare`.
The manifest includes source paths, actual/target sizes, seeds, normalization,
generation timings/statuses, and SHA256 hashes. Matching successful inputs are
reused; `--force` regenerates modules. Generation has a separate 1200-second
deadline per presentation and removes failed partial SCCs.

Use `--families torus_gaussian`, `--sizes 20 32`, or `--homology 3` to select
other inputs. `--seed 42` remains an alias for `--seeds 42`.

## Measurement protocol

- Runtime: `pruning INPUT --all-optimizations --profile --no-output`.
- Optimizations: the all-off baseline, then each of fuse I preimage, defer K
  image minimization, projected preimages, dense kernel pivots, throttled
  progress, composition generators, amortized reductions, Hom classes, phase
  quotient, and target-reduction reuse enabled alone. No switches are silently
  enabled in the baseline. Composition uses the Hom quotient internally, even
  with the separate Hom-class switch off.
- Both experiments cover every selected size and epsilon. Compact storage is
  held fixed across all configurations.
- Primary time is the CLI's internal total pruning profile, including Hom
  construction, map selection, iterations, progress formatting/writes and final
  presentation/minimization. Input loading, process startup, profile printing
  and SCC output are outside that timer. Complete process wall time is also
  recorded. Pruned SCC output is disabled.
- Epsilon is computed separately for each input and percentage from the largest
  coordinate extent over all generator and relation degrees. Different epsilon
  values are separate observations and table rows, never pooled into one median.
- Runs are sequential; configuration order is shuffled with a fixed seed.
  Synthetic seeds remain separate rows.
- Three timing observations are requested. If the first pruning time exceeds
  five seconds, only that observation is used. Failed configurations are
  attempted once. Medians, min/max, completed/attempted counts and partial
  repeats are recorded.
- Every process has a **600-second deadline**, including loading. A peak-RSS
  watchdog samples every 200 ms and stops above **12288 MiB**. This is a sampled
  guard, not a strict memory quota. Timeouts, memory stops, SIGKILL, generation
  failures and algorithm errors are distinct statuses, never fabricated timings.
- The module-pruning Hom constructor bypasses the old 100000-element conversion
  cap; the benchmark's resource limits apply instead. No alternative
  optimization configuration is substituted on failure.

## Output files

- `table.csv`: runtime or one-at-a-time optimization matrix in seconds, with
  point counts, epsilon percentages and scalar epsilons.
- `table.tex`: tables grouped by family, homology and epsilon percentage for the
  paper. Load `longtable` and `pdflscape`. Seeds have separate rows.
- `summary.csv`: one row per input/epsilon/configuration, with medians, min/max,
  sample counts and statuses.
- `observations.jsonl`: every completed attempt, flushed immediately, with the
  command, input/binary hashes, timings, iteration counts, status and limits.
- `run.json` and `build-cache-compact.txt`: experiment settings, executable hash,
  machine/Python/platform information and CMake configuration.

## Validation

After building `build/pruning`, run the small infrastructure check:

```sh
python3 -m experiments.support.check_paper_benchmarks
```

It checks no-argument entry points, all-off and individual-switch commands,
multiple epsilons on small synthetic/real inputs, repeats, automatic resume,
historical-result preservation, tables, timeouts, memory stops and LaTeX
compilation when available. It does not run the full paper experiment.
The separate `pruning_optimizations` CTest checks all 1024 small-fixture switch
combinations and exact output isomorphism.

## Earlier results

The run started on 11 November 2026 is under `paper_runs/2026-11-11/`.
Its data contains 104 inputs: 24 per synthetic family and eight observed clouds.
Those measurements and all earlier comparison protocols remain historical;
they are not relabeled or merged into the current experiment.
