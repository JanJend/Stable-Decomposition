# Research experiments

Run either experiment from the repository root, with **no arguments**:

```sh
python3 experiments/benchmark_runtime.py
python3 experiments/benchmark_optimizations.py
```

Both experiments read **`data/manifest.csv`** and reuse the same point sets and
presentations. If the shared manifest is absent, either script runs
`generate_data.py` once to prepare both size ranges. The runtime experiment uses
20, 40, 80, 160, 320, 640, 1280 and 2560 points, plus both complete large hypoxic
regions. The optimization comparison selects only Gaussian-noise tori,
Gaussian clusters and immune inputs at **80, 160 and 320 points**.
This restriction is defined in `benchmark_optimizations.py`; generation and
the shared point-set/presentation collection keep all families and sizes.
For observed clouds, each experiment selects distinct complete samples nearest
its target sizes; large hypoxic regions are added to the runtime experiment
without adding huge synthetic inputs or baseline optimization runs.

The runtime experiment enables all optimizations. The optimization experiment
compares the quotient-based baseline with each optional improvement separately
and with a full-End control. Dense and Hom algorithm comparisons are excluded.
The CLI and every benchmark configuration use A-1/2 (`optimised`) by default.
Runtime sweeps **1%, 2%, 4%, 8%** of each presentation's largest degree extent;
the optimization comparison uses only **1% and 8%**.
Settings are at the top of the benchmark scripts. Run the scripts one at a time
so their timings do not compete.

Finish data generation before starting a benchmark. The generator and benchmark
lock the shared inputs so regeneration cannot overwrite presentations during a
run. Every selected presentation is checked against the manifest before timing
starts. If inputs have changed since a previous run, rerunning creates a new
results directory and preserves the previous observations.

Results go to `results/runtime/` and `results/optimizations/`. Each contains
`table.csv` and `table.tex`, updated after every measurement; the optimization
experiment also writes `speedups.csv` and `speedups.tex`. Rerunning resumes a
matching interrupted experiment or regenerates completed tables. Changed code,
inputs or settings use a timestamped result directory, preserving earlier
measurements. Ctrl-C keeps completed observations.

## Generator-budget sweep

The optimization experiment also compares both composition methods together at
**3, 6, and 12 orbit attempts**, using hybrid seeds. Other optional optimizations
are off for these comparisons, matching the quotient-based baseline. The
one-at-a-time columns continue to use the production defaults.

```sh
python3 experiments/benchmark_optimizations.py \
  --composition-orbit-rounds 3 6 12 24 \
  --composition-orbit-seed-mode sum random hybrid \
  --composition-orbit-min-gain 1 --composition-orbit-patience 3
```

The round and seed-mode lists form a Cartesian product. `sum` uses all remaining
nonpivot maps, `random` includes each independently with probability one half,
and `hybrid` alternates these, starting with a sum. The same fixed random seed
(default 923407, set with `--composition-orbit-random-seed`) is reused across
repetitions. Patience counts consecutive attempts with zero savings or savings
below the minimum; zero patience disables early stopping. Every count is a
maximum, so runs can stop earlier. Zero rounds disables the orbit stage.

The default `--composition-orbit-pair-visits 0` automatically scales the total
allowance with the number of rounds and the End dimension. A positive value
imposes a separate cap. `--composition-sparse-seeds`,
`--composition-nilpotent-rounds`, and the remaining search budgets are also
configurable (`--help` lists all of them). Values other than the round/mode lists
apply to every sweep configuration.

`generators.csv` and `generators.tex` report selected I/K map counts, full pruning
times and baseline-relative speedups for each configuration. Complete settings,
commands and failures are saved with the observations. Changing any sweep
setting creates a new results directory; rerunning the same setting resumes it.

## Hom-size experiment

```sh
cmake --build build/local-release --target hom_reduction_benchmark
python3 experiments/benchmark_end_sizes.py
```

It uses `data/presentations/torus_gauss/torus_gauss_640_1_H1.scc` (640 points,
sample 1, RNG seed 42) with epsilon values **1%, 2%, 4%, 8%** of the largest
degree extent over generators and relations. Pass `--percentages 0.5 1 2 4`
to choose another sweep. Each epsilon runs in a separate process with a
10-minute (600-second) deadline and 12288 MiB sampled peak-RSS limit; both are
configurable. This experiment measures the global map families; it does not run
the adaptive I/K pruning iterations.

`results/End_sizes/results.csv` has **27 columns**: pointset and epsilon; serialized generator settings;
full Hom and quotient dimensions; I/K counts for **Sparse-nilpotent**,
**Complement-orbit**, and **Combined** in both spaces; the three basis/quotient
times; two shared coordinate-setup times; and six method times. Each method
time combines I and K and excludes the separately reported coordinate setup.
Combined runs Sparse-nilpotent first, then extends its certificate with
Complement-orbit; it never increases the sparse method's output count.

Both methods are enabled by default. Use `--no-sparse-nilpotent-generators` or
`--no-complement-orbit-generators` to measure only the other method. Disabled
method and Combined cells say `disabled`; they are not reported as zero time or
zero generators. The same switches work on `hom_reduction_benchmark` and the
main pruning CLI. `--no-composition-generators` remains the pruning CLI's master
switch. The one-at-a-time optimization benchmark also has separate Sparse and
Orbit columns. All `--composition-*` search controls are forwarded to the worker,
and the CSV and JSON record their values.

Failed runs retain completed measurements and mark unfinished cells `timeout`,
`memory_limit`, or `error`. Each run appends CSV rows and saves stdout/stderr
logs, detailed per-epsilon JSON records (including the selected methods), and
settings in a timestamped subfolder. An older CSV with a different layout is
preserved there as `previous_results.csv`. Run
`python3 -m experiments.support.check_end_sizes` for a small end-to-end check.

## Directory layout

```text
experiments/
  generate_data.py            Prepare shared inputs for both experiments
  plot_pointsets.py           Plot clouds and input Hilbert functions once
  benchmark_runtime.py        Runtime experiment
  benchmark_optimizations.py  Individual optimization comparison
  data/
    manifest.csv              Shared metadata; paths relative to data/
    point-sets/<family>/      One file per cloud, shared across homology degrees
    presentations/<family>/  Minimal presentations, with Hk suffixes
    pictures/                 Shared PNGs, index.csv and index.html
  results/
    runtime/                  Runtime tables and observations
    optimizations/            Optimization tables and observations
  support/                    Shared implementation and optional checks
```

## Build and prerequisites

Use Python 3.9 or newer, NumPy and Matplotlib for data generation and plotting.
Timing and reporting use the standard library. Build the pruning executable and
Hilbert renderer separately; the benchmark scripts never compile code:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGRLINA_USE_COMPACT_COLUMNS=ON
cmake --build build --target pruning experiment_hilbert --config Release -j 2
```

On Apple Silicon, add `-DCMAKE_PREFIX_PATH=/opt/homebrew` if needed. Sibling
Persistence-Algebra and AIDA checkouts supply the normal project dependencies.
Build function-Delaunay in its sibling checkout for module generation:

```sh
cmake -S ../function_delaunay -B ../function_delaunay/build -DCMAKE_BUILD_TYPE=Release
cmake --build ../function_delaunay/build --target main -j 2
```

Observed immune-cell inputs come from the sibling
[SpatialPatterningOfImmuneCells](https://github.com/MultiparameterTDAHistology/SpatialPatterningOfImmuneCells)
checkout: its `1.5mmRegions/` folder and root-level `LargeHypoxicRegion.csv` and
`LargeHypoxicRegion2.csv`. Pass `--immune-repo PATH` for another location.

## Prepare the datasets

```sh
python3 experiments/generate_data.py
```

Every synthetic family uses the same default doubling sequence starting at 20:
**20, 40, 80, 160, 320, 640, 1280, 2560**.
Every synthetic family/size/seed is a deterministic sample. Default RNG seeds
42, 43 and 44 are labelled samples 1, 2 and 3. For example, `S2_80_1.txt` has
presentation `S2_80_1_H2.scc`. Manifests record `seed` and `rng_seed` separately;
`--seeds` specifies actual RNG seeds. Sample numbers remain stable when
regenerating a subset. `--seed 42` remains an alias for `--seeds 42`.

| Family | Sampling | Default homology |
| --- | --- | --- |
| `torus_gauss` | Area-uniform torus in R3, major radius 2, minor radius 1; Gaussian sigma 0.27 and 3% clipped Cauchy outliers | H1 |
| `torus_uni` | Same mixed noise, plus replacement of 10% of points by uniform noise in [-3,3]^2 × [-1,1] | H1 |
| `gauss_clust` | Three Gaussian clusters in R3; independently jittered centers, weights 0.4/0.35/0.25, standard deviations 0.6/0.8/0.5, plus 5% clipped Cauchy outliers | H0 |
| `S2` | Uniform unit 2-sphere in R3, with the same Gaussian noise mixture and 5% clipped Cauchy outliers as clusters | H2 |
| `immuno` | Complete observed small clouds plus both complete large hypoxic regions | H1 |

The torus parameters and noise follow `../mpm_generation/noisy_torus.py`:
Gaussian perturbations on all points and random-direction Cauchy displacements
on outliers, with scale 0.6 and distances clipped to ±1.8. At least one outlier
is used in small samples. Clusters use the same displacement model with scale 2
and distances clipped to ±6; their per-coordinate variances are four times the
previous values. The `S2` family replaces the former S3 in R4.

Cluster centers are proposed with independent Gaussian jitter of standard deviation
0.75 per coordinate. Proposals are rejected unless every pair of centers is at
least 3 × (sigma_i + sigma_j) apart, keeping their three-sigma core balls disjoint.
Locations vary with the RNG seed and point count while remaining reproducible;
the three within-cluster standard deviations stay distinct. Gaussian tails and
explicit outliers may still extend beyond those core balls.

Clusters and S2 share the same noise code: Gaussian standard deviations
0.6/0.8/0.5 with weights 0.4/0.35/0.25, followed by 5% Cauchy outlier shifts
with scale 2 and distances clipped to ±6. On S2, each point independently selects
a Gaussian component; on clusters, the component also determines its center.

The cluster and S2 outlier count is max(1, floor(0.05n)): 1, 2, 4, 8, 16, 32, 64 and
128 for the default sizes. These are points selected for Cauchy displacements;
the displacements can be small, so not every selected point appears isolated.

Small real samples are selected from all cell-type subfolders by relative size
difference. Both large CSVs combine **all CD8, CD68 and FoxP3 cells**; no cell-type
filter or subsampling is applied. They contain 15,724 and 21,124 points,
respectively. CSV headers are excluded and coincident coordinates are removed.
Actual counts and source filenames are retained in the manifest.

Every family uses the **sublevel Delaunay–Čech bifiltration**: minimum enclosing
ball radius (not squared radius), and maximum vertex codensity on each simplex.
Coordinates are translated and isotropically scaled by their largest range L.
The vertex function is negative Gaussian KDE, min–max normalized to [-1,0], so
higher-density points enter first:

```text
q = (points - coordinate_minima) / L
D_i = (1/n) * sum_j exp(-||q_i-q_j||^2 / (2*h^2))
f_i = -(D_i-min(D)) / (max(D)-min(D))
```

Bandwidths are selected or specified in explicit spatial units:

- Tori: **selected separately for each cloud by leave-one-out Gaussian likelihood**,
  from 25 logarithmically spaced widths between 0.1 and 2.0 in original units.
  Larger samples generally select narrower kernels. The selected width is divided
  by L for the scaled cloud. Tori use H1 only.
- Clusters and S2: **0.2 in original units** (`add_kde.py`), so h = 0.2/L.
- Immune cells: **0.05 in scaled coordinates** (`w_multipers.py`). This uses
  the common isotropic spatial normalization here, rather than independently
  scaling x and y as that older script does.

For torus selection, each point is scored using the other n−1 points and the
proper Gaussian normalization, including the bandwidth-dependent factor h^(-3).
The mean log density selects the width; the final codensity still uses raw KDE
min–max normalized to [-1,0]. The manifest records `kde_bandwidth_method` and
`kde_bandwidth_loo_score`, alongside the actual selected scaled width.
The search is exact and blockwise, with O(25n²) time and O(128n) memory.
Scaling both distances and bandwidth by 1/L preserves Gaussian kernel weights.

`--bandwidth VALUE` overrides automatic selection or the fixed family value
in the units above. The final KDE includes self-contributions and outliers;
the common Gaussian normalization factor cancels under min–max normalization.
A constant D gives f=0. KDE is exact and blockwise, using O(n²) time and O(128n)
memory, including for the large observed clouds. The manifest records the
actual scaled `kde_bandwidth`; multiply by `spatial_scale` for original units.

Use `--families torus_gauss`, `--sizes 20 40`, or `--homology 1` to prepare a
subset. Immune generation always includes both large regions. Matching
successful presentations are reused; `--force` regenerates them. The generation
deadline is 1200 seconds per presentation, with partial failed SCCs removed.
Source paths, sizes, seeds, normalization, generation statuses/timings and SHA256
hashes are recorded in the shared manifest. `--output-dir` supports isolated
checks without creating another experiment-specific data store.

## Inspect the input modules

```sh
python3 experiments/plot_pointsets.py
```

Open [`data/pictures/index.html`](data/pictures/index.html) for the shared gallery
of point sets and Hilbert functions before pruning. Images live under family
subfolders: `NAME_points.png` for a cloud and `NAME_Hk.png` for its Hilbert
function. A point set used in H1 and H2 is plotted once. Matching images are
reused. `index.csv` records plot bounds, bandwidths and input hashes.

Immune clouds are plotted in 2D. Clusters, tori and S2 use **3D plots**, with
normalized KDE as colour and every point included. Use `--experiment runtime`
or `--experiment optimizations` to render a subset into the same directory;
the default gallery includes all shared inputs. `--data-dir PATH` supports an
isolated generated bundle; `--hilbert-renderer PATH` selects another built
renderer. `--self-check` exercises the 2D/3D/legacy-4D point renderer.

Hilbert images show dim H_k(r,t) over F2 at 1000 × 700 coordinates, using the
complete projective resolution, including syzygies. X is radius and Y is
negative normalized KDE. White is zero; positive dimensions use logarithmic
blue colours. Each image has its own radius extent and colour maximum. Very
narrow features may fall between sample coordinates. Images belong to the
input data and are shared across all epsilon values and optimization settings.

## Measurement protocol

- Runtime: `pruning INPUT --profile --no-output` uses all optimizations and Hom Algorithm A-1/2 by default.
- Optimizations: Base uses the Hom quotient and Algorithm A-1/2 with all five
  optional improvements explicitly disabled. Proj, Amort, Comp, Adaptive and
  Reuse each enable only their named improvement. The Full End control instead
  passes `--full-end-space`, computes the full End_(2 epsilon) basis without
  quotienting by shifted End_0, and visits that basis in the I/K loops.
  Its absolute time is labelled Full End; the percentage column is labelled
  Quot and reports `100*(Full End/Base - 1)`, measuring how much faster the
  quotient-based standard is. Other percentages use `100*(Base/Variant - 1)`.
  Both timings include basis preparation and all pruning iterations.
  Dense pivots, fused I preimages, deferred K image minimization and
  percentage-based progress are unconditional and have no comparison columns.
  Hom algorithm comparisons are excluded; every benchmark configuration uses
  A-1/2. Use `--hom-algorithm=hom-exactness` outside the benchmark to select B.
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
  Individual seeds remain in `summary.csv` and `observations.jsonl`. Optimization
  tables average the three per-seed median runtimes with equal weight; immune
  rows use their single observed sample. Speedups are ratios of these mean
  runtimes, not averages of percentages. A missing or failed seed never silently
  drops out of a mean. Runtime tables retain separate rows for each seed.
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

- `table.csv`: runtime matrix in seconds, or optimization mean runtimes with
  one row per family, point count and epsilon percentage. Per-seed timing
  medians receive equal weight; unequal repeat counts do not weight a seed more.
- `table.tex`: optimization results use one compact, standard `table`/`tabular`
  with point-set, size and epsilon columns, followed by the configuration times.
  There are no seed or homology columns, landscape pages, page-size changes or
  resizing. The complete optimization selection gives 18 rows. Runtime keeps
  its detailed per-input tables and requires `longtable`.
- `speedups.csv` and `speedups.tex`: optimization speedups from the mean runtimes.
  Quot reports `100*(Full End/Base - 1)`; other columns report
  `100*(Base/Variant - 1)`. A 2x speedup is +100%; a runtime twice the reference
  is -50%. The zero Base column is omitted. The LaTeX output uses the same
  compact layout and 18 averaged rows. CSVs retain full precision; LaTeX cells
  show three significant digits. Pending cells are `--`, failed comparisons
  are `n/a`, and `*` marks partial timing repeats. Per-seed failures and counts
  remain in `summary.csv`. Runtime has no percentage comparison.
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

After building `experiment_hilbert` and function-Delaunay, check the data and
plotting pipeline with:

```sh
python3 -m experiments.support.check_data
```

This checks 3D samples, mixed outliers, exact KDE, automatic torus bandwidth
selection and overrides, all immune cell types, shared input selection, H2,
presentation reuse and picture output.

It checks no-argument entry points, explicit baseline and individual-option commands,
multiple epsilons on small synthetic/real inputs, repeats, automatic resume,
historical-result preservation, tables, timeouts, memory stops and LaTeX
compilation when available. It does not run the full paper experiment.
The separate `pruning_optimizations` CTest checks all 128 small-fixture switch
combinations and exact output isomorphism.
