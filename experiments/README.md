# Experiments

## Run

Run commands from the repository root. Requirements:

- Python 3.9+ and NumPy (`python3 -m pip install numpy`).
- CMake 3.28+, a C++17 compiler and Boost (timer, graph, chrono).
- Sibling checkouts: `../Persistence-Algebra`, `../AIDA` and
  `../function_delaunay` (with CGAL, Eigen and GUDHI installed).
- For immune inputs, the sibling
  [SpatialPatterningOfImmuneCells](https://github.com/MultiparameterTDAHistology/SpatialPatterningOfImmuneCells)
  checkout containing `1.5mmRegions/`, `LargeHypoxicRegion.csv` and
  `LargeHypoxicRegion2.csv`.

Build once:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGRLINA_USE_COMPACT_COLUMNS=ON
cmake --build build --target pruning hom_reduction_benchmark -j 2
cmake -S ../function_delaunay -B ../function_delaunay/build -DCMAKE_BUILD_TYPE=Release
cmake --build ../function_delaunay/build --target main -j 2
```

On Apple Silicon, add `-DCMAKE_PREFIX_PATH=/opt/homebrew` to the configure
commands if dependencies are not found.

Prepare the shared inputs, then run the experiments one at a time:

```sh
python3 experiments/generate_data.py
python3 experiments/benchmark_runtime.py
python3 experiments/benchmark_optimizations.py
python3 experiments/benchmark_end_sizes.py --worker build/hom_reduction_benchmark
python3 experiments/end_quotient_size.py --worker build/hom_reduction_benchmark
```

To profile only the **1280-point Gaussian torus H1** (sample 1, RNG seed 42),
using the runtime benchmark's all-optimizations path with larger torus search budgets:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DGRLINA_USE_COMPACT_COLUMNS=ON
cmake --build build --target pruning_torus_runtime -j 2
python3 experiments/benchmark_torus_runtime.py
```

This uses the production pruning entry point, including fast I, the direct
single-solve quotient, compact columns and automatic worker detection. It
reuses the manifest's presentation without generating inputs. Defaults are
1%, 2%, 4% and 8% epsilon, one fresh process per epsilon, an 1800-second deadline
and a 12288 MiB memory limit. For one measurement:

```sh
python3 experiments/benchmark_torus_runtime.py --percentages 1 --timeout 1800
```

Use `--threads 1` for serial measurements or `--threads N` for a chosen worker
budget; `--repeats N` retains each repetition separately. Add `--standard-i`
to compare the fixed-per-pass I rule with fast I. Add `--parallel-i` to fix J
per pass and process groups of homomorphisms concurrently, then intersect their
results. With a 12-thread budget, I and K each receive 6 workers. A one-thread
budget falls back to the serial fixed-target loop. J is formed without a separate
generator minimization; `--minimize-i-target` restores that work for comparison.
For separate archives of the three current variants:

```sh
python3 experiments/benchmark_torus_runtime.py --percentages 1 --threads 12 --output-dir experiments/results/torus_i_fast
python3 experiments/benchmark_torus_runtime.py --percentages 1 --threads 12 --standard-i --output-dir experiments/results/torus_i_fixed_serial
python3 experiments/benchmark_torus_runtime.py --percentages 1 --threads 12 --parallel-i --output-dir experiments/results/torus_i_fixed_parallel
```

The composition-search options
from the other experiments are supported too. This torus experiment defaults
to 262144 certificate visits and a 262144-unit nilpotent-reduction budget per
side, shared by certificates and products. These are experiment-specific
overrides; ordinary pruning uses 8192 for both options. To use the ordinary
production budgets in this torus experiment:

```sh
python3 experiments/benchmark_torus_runtime.py --percentages 1 --composition-certificate-visits 8192 --composition-nilpotent-products 8192 --output-dir experiments/results/torus_runtime_1280_original_budget
```

The production path stores each map as one sparse coefficient column and lifts
only the active map. Restriction selection processes batches of 32 maps using
shared target reducers and one incremental independence reduction. Map order,
selected indices and the algebraic stopping criteria are preserved. Fast I
reuses its image and prepared target while I's embedded generators remain
identical; a different representation conservatively triggers a refresh.

Results are under
`experiments/results/torus_runtime_1280/`:

- `results.csv`: status, internal pruning time, process wall time, peak RSS,
  map counts, iterations and output sizes for each epsilon/repetition.
- `steps.csv`: every timed operation call, with start/end times, duration and
  I/K pass number. Unfinished operations have blank durations. Pass 0 denotes
  setup, finalization or an aggregate phase timer.
- `timings.csv`: completed calls and summed durations by operation, plus
  unfinished-call counts.
- A timestamped archive contains these CSVs, settings, input/binary hashes,
  build configuration, the completed output as `.output.scc`, full stderr and raw `.trace.tsv` start/end events,
  flushed during execution. Completed steps survive a timeout or interruption.

Setup dimensions, map counts and the peak RSS reached by setup are emitted
before I/K begin, so these remain available if a later phase fails. On a failed
run, recorded peak RSS is a lower bound on the process's eventual peak.

The setup breakdown separates End_0, the direct Hom quotient and composition
generator selection (including its coordinate setup). `map_lifts` records
coefficient-family setup; temporary map lifting is included in I/K phases. I/K
operations include restriction selection, target preparation, images,
preimages, reductions and convergence checks.
`i_shift` measures forming J without minimization; `i_target_reduce` measures
the optional separate minimization. `i_intersection` measures merging the
parallel workers' submodules. Parallel I operation totals sum durations across
workers and can exceed the I phase's wall time. Loading/sorting is timed
separately; pruning excludes input loading and SCC output, as in the runtime
benchmark. Timings include trace overhead. `total` contains the other pruning
stages, `basis` contains the End/quotient/composition breakdown, and I/K phases overlap when
multiple workers are used; do not add nested or overlapping rows.

Run the small regression check with
`python3 -m experiments.support.check_torus_runtime`.

To measure how many generators **I and K hold throughout pruning**, relative to
the presentation of X:

```sh
cmake --build build --target pruning_torus_runtime -j 2
python3 experiments/benchmark_generator_sizes.py
```

This reuses the production worker and the shared sample-1 Gaussian torus H1,
Gaussian clusters H0 and complete immune H1 inputs at **640 and 1280 target
points**. Actual observed point counts are retained (the 1280 immune input has
1281 cells). Defaults are **0.5%, 1% and 2%** epsilon, the standard I rule,
automatic workers, a 600-second deadline and a 12288 MiB memory limit per run.
Use `--fast-i` to measure the fast I rule instead. Other examples:

```sh
python3 experiments/benchmark_generator_sizes.py --sizes 640 1280 2560 --percentages 1 2
python3 experiments/benchmark_generator_sizes.py --families torus_gauss --sizes 1280 --timeout 1800
python3 experiments/benchmark_generator_sizes.py --input test_presentations/test8.scc --percentages 10 --threads 1
```

Results are saved under `experiments/results/generator_sizes/`:

- `results.csv`: input generators and relations, peak stored and peak reduced
  I/K generator counts, final counts, iterations, runtime, status and ratios.
  `i_peak_generators`/`k_peak_generators` include unreduced candidates;
  `i_peak_reduced_generators`/`k_peak_reduced_generators` exclude them.
- `generators.csv`: initial counts, each I update and K preimage before/after
  minimization, pass endings and final counts. Separate J and L rows track
  `J=eta(I)` and `L=eta(K)+sum f(K)`, including L's unreduced accumulation batches.
  `map` is the 1-based index in the phase's full map family, before adaptive
  restriction; 0 means no individual map. Pass 0 denotes initialization or the Q=0 shortcut.
- Each count is divided by both `generators(X)` and
  `presentation_size=generators(X)+relations(X)`. Zero denominators give blank
  ratios. Counts are columns of the actual generating families, with no extra
  minimizations or submodule presentations. The initial I uses the supplied
  generators of X; internal kernel workspaces are not measured.
- A timestamped archive retains settings, hashes, raw flushed traces, stderr
  and CSVs. Timeouts/interruption preserve completed events; their peaks are
  labelled lower bounds and missing counts stay blank. The algorithm, including
  final quotient formation, runs in a fresh process for every input/epsilon.

Run the small regression check with
`python3 -m experiments.support.check_generator_sizes`.

For the **0.5% generator-growth curves through 640 points**, rerun all six
target sizes with a **one-hour deadline per input**:

```sh
python3 experiments/benchmark_generator_sizes.py --sizes 20 40 80 160 320 640 --percentages 0.5 --timeout 3600 --output-dir experiments/results/generator_sizes_to640_eps0p5_1h
python experiments/plot_generator_sizes.py
```

Use a Python with Matplotlib for the second command. It saves a PNG, an SVG
and `percentages.csv` beside the results. The three panels each show two curves:
the peak reduced generator counts of I and K, divided by
`generators(X)+relations(X)` and multiplied by 100. Both axes are logarithmic,
and the x-axis uses actual point counts (including the 41-cell immune sample
near target 40). The panels share their percentage scale; the dotted line
marks 100%. Unfinished phases appear as lower-bound arrows, and zero values
are explicitly noted rather than replaced by artificial positive values.
Supply `--input PATH/results.csv` to plot an archived run. Check the plot with
`python -m experiments.support.check_generator_size_plot`.

To run the End-size benchmark only on one immune sample:

```sh
python3 experiments/benchmark_end_sizes.py --immuno-only 640
python3 experiments/benchmark_end_sizes.py --immuno-only 1280
```

These select `T_E_ROI_40_locations_CD8` (640 cells) and
`T_A_ROI_36_locations_CD8` (1281 cells), respectively, both in H1. They use
0.5%, 1% and 2% epsilon and save counts and timings separately under
`experiments/results/End_sizes_immuno_640/` or `End_sizes_immuno_1280/`.
The existing `--immuno-640-only` spelling remains an alias for `--immuno-only 640`.

For the three predetermined 640-point examples (Gaussian torus H1, Gaussian
clusters H0 and immune H1), compute dimensions and left/right E-generators of Q
without running the later pruning iterations:

```sh
python3 experiments/benchmark_end_sizes.py --640-only --generators-only
```

This writes `dimensions.tex`, `generators.tex` and their combined `results.tex`
under `experiments/results/End_sizes_640/`, with CSV data, settings and archived
logs. E is End_0; both generating families generate Q under its respective
E-action. Their cardinalities are certified upper bounds, not minimum generator
numbers. To select one or two of these inputs, add, for example,
`--families torus_gauss gauss_clust` or `--families immuno`. The existing
`--input PATH.scc` option selects one arbitrary presentation.

To time every restriction reduction on the 640-point Gaussian torus, including
the later pruning passes:

```sh
python3 experiments/benchmark_end_sizes.py --640-only --families torus_gauss --time-steps --output-dir experiments/results/End_sizes_torus_640_timed
```

`--time-steps` defaults to off. When enabled, `restrictions.csv` has an additional
`seconds` column for each I/K pass's restriction selection, including zero-retained
passes. These times measure restriction and independence reduction, excluding
the subsequent submodule updates. `timings.csv` already records the main stages
(End spaces, quotient, coordinate setup, combined left/right generator selection),
total restriction-selection times, and full I/K phase times. Restriction-selection
times are included in the phase times; do not add them again. The worker accepts
the same `--time-steps` flag. Existing stage timers remain enabled in the benchmark.

The End-size benchmark uses the explicit quotient route. Its `Quotient seconds`
column measures normalization of the shifted End_0 image and quotient reduction;
it excludes the separately timed End_0 and End_(2 epsilon) computations. For a
comparison with the direct single-solve Q routine on the same 640-point torus,
run the existing paired benchmark:

```sh
cmake --build build/local-release --target end_quotient_benchmark -j 2
python3 experiments/benchmark_end_quotient.py --sizes 640 --families torus_gauss --repeats 1 --output-dir experiments/results/end_quotient_torus_640_comparison
```

This records both methods in `results.csv`, using a fresh process for each.
`explicit_seconds` includes End_0, shifted End and quotient formation;
`explicit_with_E_seconds` excludes End_0, which is useful when E is already
available for module actions. `end_d_0_seconds` measures the direct Q computation
without an End_0 basis. Dimensions must agree before speedups are reported.

For curves on the three **320-point** examples, sweep epsilon from **0% to 5%**
inclusive in **0.1%** steps (51 measurements per input):

```sh
python experiments/end_quotient_curve.py
```

Use a Python environment with Matplotlib installed (`python` in the local
Anaconda environment already has it). This saves all 153 measurements to
`experiments/results/end_quotient_curve/results.csv` and opens one figure with
Gaussian torus, clusters and immuno panels. Each panel shows
`dim End_(2 epsilon)` and `dim Q`, with `Q=End_(2 epsilon)/(eta End_0)`.
Save the figure manually using the Matplotlib toolbar; no plot file is written.
The inputs and per-axis percentage convention match the End quotient size
experiment. `End_0_dimension` and actual epsilon coordinates are also in the CSV;
settings, raw observations and worker logs are archived per run. The CSV is
updated after every measurement; failed or pending dimensions stay blank and
appear as gaps in the plot. To compute without a window or reopen the CSV later:

```sh
python experiments/end_quotient_curve.py --no-plot
python experiments/end_quotient_curve.py --plot-only
```

For the same sweep on the three **640-point** examples:

```sh
python experiments/end_quotient_curve.py --points 640
python experiments/end_quotient_curve.py --points 640 --plot-only
```

The 640-point CSV and archived logs are saved separately under
`experiments/results/end_quotient_curve_640/`. The figure uses the same titles
and legend, with no y-axis label.

Compare the explicit quotient computation with the new single-solve `End_d_0`
routine on the precomputed Gaussian torus H1, clusters H0 and immune H1 inputs:

```sh
cmake -S . -B build/local-release -DCMAKE_BUILD_TYPE=Release
cmake --build build/local-release --target end_quotient_benchmark -j 2
python3 experiments/benchmark_end_quotient.py
```

The default sweep doubles target sizes from 20 through 2560 and uses 0.5%, 1%
and 2% epsilon. For a shorter comparison, add `--sizes 80 160 320`,
`--families torus_gauss immuno`, or `--percentages 1`. Inputs are synthetic
sample 1 (RNG seed 42) and complete observed clouds; actual point counts are
recorded separately. No inputs are regenerated.

Each method runs in a fresh process with A-1/2. The explicit route computes
End_0 and then calls the supplied-E overload; the new route calls
`End_d_0_coefficients(eta, algorithm)` without an E basis. Internal times exclude
loading, sorting, initial row caches and shift setup. `speedup_explicit_over_new`
is explicit total time divided by new time; values above one favor the new
routine. `speedup_with_E_over_new` compares the explicit time after End_0 is
available against new time. The latter is useful when E is already needed for
composition. Both ratios are left blank on a failure or dimension mismatch.

Results are saved after each completed pair to
`experiments/results/end_quotient_comparison/results.csv`; timestamps archive
settings, hashes, raw repetitions and process logs. The script checks shifted
Hom, quotient and shift-image dimensions. It reports medians of three paired
runs, alternating method order; a first pair with either timing over five seconds
is used alone. Defaults are 600 seconds and 12288 MiB per process, configurable
with `--timeout` and `--memory-limit-mb`. Interrupting preserves completed pairs.
Run the small regression check with
`python3 -m experiments.support.check_end_quotient_benchmark`.

Generation reuses matching presentations. Pass `--immune-repo PATH` to the
generator if the immune checkout is elsewhere. Benchmark sizes, epsilon values,
repeats and limits are defined at the top of the benchmark scripts; `--help`
lists available command-line overrides.

Results are under `experiments/results/`:

- `runtime/`: `table.csv` and `table.tex`.
- `optimizations/`: `table.*`, `speedups.*` and generator-budget `generators.*`.
- `End_sizes/`: `results.csv` and `results.tex`, plus separate `timings.csv` and
  per-pass `restrictions.csv`, with an archive for each run.
- `end_quotient_size/`: `results_320.csv`, `results_640.csv`, `results_1280.csv`
  and matching `.tex` files; `results.tex` contains all three tables. Each has
  three dimension rows and nine data columns (torus, clusters, immuno, each
  at 0.5%, 1% and 2%). Settings, raw observations and logs are archived per run.
- `end_quotient_curve/`: `results.csv`, with 51 epsilon values per 320-point
  input, plus an archive for each run. Plots open interactively and are not saved.
- `end_quotient_curve_640/`: the same sweep and CSV layout for 640-point inputs.

Raw observations and settings are saved alongside the tables. Runtime and
optimization runs resume matching interrupted runs; changed inputs, binaries or
settings preserve earlier results in separate directories. The size experiment
starts a new archived run on every invocation.

## What is measured

All experiments use Hom Algorithm A-1/2. Epsilon scales each bounding-box side over generator and relation degrees
separately: `epsilon = (P/100 * Lx, P/100 * Ly)`. Raw epsilon and shift values
are vectors. Results from the previous largest-side convention remain in their
existing files; rerun an experiment to measure the new convention.

- **Runtime:** pruning with all optimizations, over synthetic point counts
  20, 40, 80, 160, 320, 640, 1280 and 2560, plus both complete large hypoxic
  regions; epsilon values are 1%, 2%, 4% and 8%.
- **Optimizations:** the quotient-based baseline versus each optional improvement
  separately and a Full End control. The **Direct Q** column measures the
  single-solve `End_d_0` routine against a baseline computing End_0 and the
  shifted Hom separately. The baseline and other individual variants use
  `--no-single-solve-quotient`; normal runs and the all-optimizations runtime
  experiment keep it enabled. This remains `single-solve-quotient` in CSV fields
  and command-line options. With composition enabled, End_0 is computed separately
  for the actions in the direct route and retained in the old route; the timed
  pruning includes that cost. **Fast I** enables `--update-i-target`:
  recompute the shift image of the current I before each map. Base and the other
  individual variants keep this image fixed for the whole pass
  (`--no-update-i-target`). Both modes use the same K loop. With adaptive
  restriction enabled, map selection still uses I at the start of each pass;
  target reductions are refreshed whenever the target changes. The runtime
  experiment enables fast I along with the other optimizations.
  Inputs are Gaussian-noise tori, Gaussian
  clusters and immune samples near 160, 320 and 640 points, at 0.5%, 1% and 2% epsilon.
  A separate generator-budget comparison combines Nilpotent reduction and
  Complement-orbit with 12 and 24 orbit attempts and hybrid seeds.
- **End sizes:** exactly seven point sets, independent of the optimization
  selection: Gaussian-noise tori, immune samples and Gaussian clusters at
  target sizes **640 and 1280**, plus complete **LargeHypoxicRegion H1**.
  Synthetic inputs use sample 1 (RNG seed 42). The 1280-target immune cloud
  has 1281 actual cells; it remains complete. Epsilon values are 0.5%, 1% and 2%,
  giving 21 rows without seed averaging. The table shows `dim End_0`,
  `dim End_(2 epsilon)`, `dim Q` for `Q=End_(2 epsilon)/(eta End_0)`, the
  selected left/right generator counts `|B_l|` and `|B_r|` after both methods,
  and the maximum retained I/K restriction counts over actual pruning passes.
  `restrictions.csv` lists every pass, including initial and terminal zeros;
  timings are separate. The worker reuses its computed quotient and selected
  families during those loops. End_0 supplies the action on Q; its basis is
  not replaced by B_l or B_r. The implementation note is
  `explanations/budgeted_composition_generators.tex` (and its compiled PDF).
  Nilpotent reduction now first certifies single-row/column quotient maps by
  elementary lowering factorizations, inserting their known unit coordinates
  without multiplication or a coordinate solve. Known valid factors are tried
  first; up to 128 additional candidates per side receive cached relation checks.
  Quotient coordinates are prepared only when a product or a newly constructed
  factor's membership check needs them. Both factors must be module maps;
  a new intermediate factor must also belong to the supplied Hom subspace.
  When direct certificates fill the nilpotent budget, the small ordinary-seed
  trials are skipped; Complement-orbit still runs. Successful certificates
  share the nilpotent-reduction product budget. Set
  `--composition-certificate-visits 0` to disable this pass for comparisons
  (default cap: 8192 candidates per side).
  With ordinary defaults, certificates and sampled lowering products share
  4096 units per side, and certificates plus all nilpotent-method products
  share 8192 units. The independent Complement-orbit allowance is unchanged.
  All three current explanatory notes include the same
  [production algorithm](../explanations/composition_generators_implementation.tex),
  covering candidate order, relation checks, rollback and stopping.
  On the fixed 1280-point torus at 1% of each bounding-box side, the tuned
  default took a median 0.110985 seconds over three runs, versus 35.986382
  seconds for the preceding certificate version, including coordinate setup.
  It kept 89238 maps on each side, versus 89176/89181 previously. Computing
  E/Q and later pruning are outside those times. Raw measurements are in
  `results/nilpotent-1280-2026-11-25/`.
  The production benchmark's `quotient_coordinates` row is now zero; any
  required coordinate preparation is included in `composition_generators`.
  The standalone worker's `--compare-methods` mode retains the older separate
  method and full-space comparisons.
- **End quotient size:** Gaussian-noise tori (H1), Gaussian clusters (H0) and
  immune samples (H1) at target sizes **320, 640 and 1280**, using synthetic
  sample 1 and the complete observed clouds. Epsilon values are 0.5%, 1% and 2%.
  The worker's `--dimensions-only` mode stops after `dim End_0`,
  `dim End_(2 epsilon)` and `dim Q`, skipping generator selection and pruning.
  Each invocation starts a new archived run and updates the three size tables
  after every measurement; completed dimensions survive a later timeout.

Inputs share `experiments/data/manifest.csv`: noisy tori (H1), Gaussian clusters
(H0), a noisy 2-sphere (H2) and observed immune clouds (H1). Synthetic RNG seeds
are 42, 43 and 44; observed clouds are used whole. Presentations use the sublevel
Delaunay–Čech bifiltration with radius and negative normalized Gaussian KDE.
Sampling and bandwidth details are in [generate_data.py](generate_data.py).

Runtime and optimization tables report internal pruning time, including Hom
construction and final minimization, excluding input loading and SCC output.
Process wall time is recorded separately. Up to three repetitions are summarized
by their median; a first run over five seconds is used alone. Optimization tables
average these medians equally across seeds; runtime keeps seeds separate.
Speedups are `100 * (Base / Variant - 1)`; Quot uses
`100 * (Full End / Base - 1)`. Each process has a 600-second deadline, except
LargeHypoxicRegion H1 in the End-size experiment: its deadline is **60 minutes
per epsilon**, configurable with `--large-timeout` (seconds). The sampled memory
limit is 12288 MiB. Failures are recorded as statuses, not timings.
