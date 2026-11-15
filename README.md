
# Multiparameter Persistence Pruning

A C++ implementation for computing the pruning of a multiparameter persistence module.

## Overview

This tool reads an SCC module, computes its pruning using the module framework,
and writes a minimal SCC presentation. Optional flags generate Hilbert function
images and compare the input/output decompositions using AIDA.

## Building
From this directory:
```bash
cmake -S . -B build
cmake --build build --target pruning -j 2
```
Compact columns are **on by default** for graded matrices. To use ordinary
`std::vector` columns, configure with `-DGRLINA_USE_COMPACT_COLUMNS=OFF`;
use `ON` to switch back, then rebuild. This is a compile-time storage choice,
applied consistently to pruning, embedded AIDA, and tests. Some ungraded scratch
matrices continue to use vectors in either build. `--help` and `--profile` report
the compiled column storage.

Use `-DCMAKE_BUILD_TYPE=Debug` for debug symbols; the executable is then named
`pruning_debug`. Build all targets and run the tests with:

```bash
cmake --build build -j 2
ctest --test-dir build --output-on-failure
```

Module fixtures live in `test_presentations/`, including the hand-checked
examples in `test_presentations/tests_with_solutions/`. The single
`tests/pruning_solutions_test.cpp` regression computes their pruning and checks
exact isomorphism against the supplied solutions:

```bash
cmake --build build --target pruning_solutions_test -j 2
ctest --test-dir build -R '^pruning_solutions$' --output-on-failure
```

Other functional tests and their instructions live in
[`tests/functional/`](tests/functional/README.md).

Synthetic and observed immune-cell runtime experiments, generation scripts and
measurement tables live in [`experiments/`](experiments/README.md).
Its [command guide](experiments/README.md) explains which
script to run; internal helpers and worker sources live in `experiments/support/`.

The default workspace layout has sibling `Persistence-Algebra` and `AIDA`
checkouts. The build reuses AIDA's bundled CLI11 parser. Persistence-Algebra
bundles the PNG writer used for Hilbert images; no Python or GUI library is needed.
Their locations can be overridden with `PERSISTENCE_ALGEBRA_DIR`, `AIDA_DIR`,
and `CLI11_INCLUDE_DIR`. Boost is required.
`-DPRUNING_WITH_AIDA=OFF` omits the AIDA library integration.

On Apple Silicon, Boost must also be built for ARM. If the linker reports
`found architecture 'x86_64', required architecture 'arm64'` for Anaconda's Boost,
select an ARM Homebrew Boost installation and clear the cached Boost paths:

```bash
cmake -S . -B build -U '*[Bb]oost*' -DBoost_ROOT="$(/opt/homebrew/bin/brew --prefix boost)"
cmake --build build -j 6
```

The `Could NOT find OpenMP` messages are harmless; OpenMP is optional.

`main.cpp` shows the workflow: read the module, choose epsilon, prune, and select
which outputs to write. Compiled implementations live in `src/` and
declarations in `include/`. `cli.hpp`/`cli.cpp` handle output, plotting, and AIDA;
`utils.cpp` defines the options. CMake shares the `stable_pruning` library between
the executable and tests. General algebra now lives in Persistence-Algebra,
with compatibility adapters here. See [the extraction review](docs/algebra-extraction.md)
for the full function mapping, behavioral notes, and tests.

## Usage

```bash
./build/pruning <file_path> [options]
```

### Options

- `-e, --epsilon <value>` - Set epsilon (default: 1% of the largest coordinate extent of all presentation degrees; fallback 0.01).
- `--delta <value>` - Compatibility alias for `--epsilon`.
- `-o, --output <path>` - Choose the output SCC path; missing directories are created.
- `--hilbert` - Save Hilbert PNGs for the input and pruned module.
- `-v, --version` - Display the version without running pruning.
- `--image-size <pixels>` - Set the heatmap grid size (default 500; range 128–2048, plus margins for labels).
- `--aida` - Run AIDA directly on both modules and save their decompositions and comparison.
- `--old` - Use the original matrix pruning instead of the module implementation.
- `--quick` - With `--old`, enable the original quick algorithm.
- `--compare` - Run matrix pruning followed by module pruning on independent copies and compare timings and output invariants. AIDA checks run automatically when built in. Both use `quick=false`; the module result is saved. Cannot be combined with `--old`.
- `--no-output` - Skip the pruned SCC; requested images and AIDA results are still saved.
- `--memory-limit-mb` - Optional sampled peak-RSS limit in MiB; zero disables it.
- `--profile` - Print module pruning times, call counts, iteration counts and enabled optimizations.
- `--no-projected-preimages` - Disable tracking only required kernel coordinates in two-parameter preimages.
- `--no-composition-generators` - Disable both composition reductions (master switch).
- `--no-sparse-nilpotent-generators` - Disable nilpotent actions and small sparse candidates.
- `--no-complement-orbit-generators` - Disable orbit compression of sums of remaining maps.
- `--no-amortize-reductions` - Reduce the K image sum after every image instead of batching up to 16 images.
- `--hom-algorithm hom-exactness|optimised|full-restriction` - Select Algorithm B, A-1/2 (default), or A for both Hom computations.
- `--full-end-space` - Use the full End_(2 epsilon) space without quotienting by shifted End_0.
- `--no-adaptive-restriction` - Disable removing dependent map restrictions in each I/K pass.
- `--no-reuse-target-reductions` - Disable reusing the reduced fixed I-pass target generators.
- `--all-optimizations` - Enable all five optional optimizations and the Hom quotient (already the default), keeping the selected Hom algorithm.
- `-h, --help` - Show usage and all options.

Options may appear before or after the input. `--epsilon=0.05` also works.

All five optional optimizations default to **on**; use the `--no-...` flags to
turn them off. The corresponding positive flags remain available. Fused I
preimages, deferred K image minimization and percentage-based progress updates
and indexed kernel pivot lookup are unconditional and have no switches.
Hom defaults to **Algorithm A-1/2**; `--hom-algorithm=hom-exactness` selects B. Optimization flags and `--profile`
apply to module pruning and cannot be combined with `--old` or `--compare`.
The legacy matrix algorithms are unchanged.

```bash
# Default optimized run, then disable projected preimages
./build/local-release/pruning input.scc -e 0.05 --no-output --profile
./build/local-release/pruning input.scc -e 0.05 --no-output --profile --no-projected-preimages

# Paper experiments: run without arguments
python3 experiments/benchmark_runtime.py
python3 experiments/benchmark_optimizations.py
```

Build `build/pruning` separately with CMake in Release mode with compact columns.
Both scripts run that executable and sweep point counts and epsilon percentages on synthetic
and observed real inputs. Runtime uses all optimizations; the optimization
comparison uses A-1/2 and the Hom quotient with the optional improvements
explicitly disabled, then enables each improvement alone. Full End omits the
quotient; the Quot percentage measures the baseline's speedup over that control.
Dense pivots, fused I, deferred K minimization and throttled progress run in every configuration.
Sparse immediately precedes Adaptive; Orbit has its own column; percentage tables omit the zero Base column.
The optimization comparison uses Gaussian tori, clusters and immune inputs at
80, 160 and 320 points, with epsilon 1% and 8%. Its compact portrait tables
average the three synthetic seed medians and omit seed and homology columns.
Runtime retains all families and the full 20-to-2560 size range. Data generation
also retains every family and size.
CSV and LaTeX tables are written automatically, and rerunning continues saved
results. Change the few experiment settings at the top of each Python file.
See [experiments/README.md](experiments/README.md) for data preparation and output paths.

The standalone [`hom_reduction_benchmark.cpp`](hom_reduction_benchmark.cpp) reports
set sizes and incremental times for `End(X)_(2 epsilon)`, its quotient by shifted
`End(X)_0`, left/right composition generators of that quotient, and left/right
generators of the full Hom space without that quotient. It accepts any bigraded
SCC module presentation and a list of epsilon values. Values use the input's
coordinate units by default; `--percent` interprets them as percentages of the
largest bounding-box side over all generator and relation degrees:

```bash
cmake --build build/local-release --target hom_reduction_benchmark
./build/local-release/hom_reduction_benchmark input.scc 0 0.005 0.01 0.02 --adaptive
./build/local-release/hom_reduction_benchmark input.scc 1 2 4 8 --percent --adaptive
```

For percentage `p`, `epsilon = p/100 * max(xmax-xmin, ymax-ymin)` and the
pruning shift is `(2*epsilon, 2*epsilon)`. A ready-to-run 640-point Gaussian-torus
H1 sweep is `python3 experiments/benchmark_end_sizes.py`; it saves a CSV and logs
under `experiments/results/End_sizes/`.

`End(X)_0` is computed once and its time is included in each route total. Stage
times exclude loading, sorting, printing and the copy needed to preserve the
full Hom basis for comparison. Coordinate setup is shared and timed separately
from Sparse-nilpotent, Complement-orbit, and Combined selection. These searches are certified; they need not be minimal. Optional `--adaptive` runs the pruning
pair and reports candidate and retained map visits summed over the actual I/K
passes, with selection times and repeated global-basis preparation time shown
separately. Tables go to stdout and library progress goes to stderr. A small
check with known dimensions can be run from the repository root:

```bash
cmake -DBENCHMARK="$PWD/build/local-release/hom_reduction_benchmark" -P tests/functional/hom_reduction_benchmark.cmake
```

The benchmark uses the module's internal total timer, excluding input preparation
and file I/O. It captures progress output identically for every variant and
shuffles run order with a fixed seed. Thus it includes progress formatting and
log-file writes, but does not measure rendering in an interactive terminal. It
reports timings, not output-equivalence checks. The `pruning_optimizations` CTest
checks all 128 combinations against the baseline: both submodules of the pruning
pair, iteration counts, and exact isomorphism of the output, including grid and
CSC boundaries. Its executable optionally accepts `FILE EPSILON` to check an
additional input using 64 representative combinations (all Hom/composition/adaptive/reuse
combinations with four kernel/batching settings). Kernel-related switches affect the I/K iteration operations;
the final presentation/minimization retain their existing paths. The composition
switch separately changes construction of the map families.
In other parameter dimensions, the kernel switches retain the existing kernel
implementation.

The C++ API accepts the same options as its last argument:

```cpp
stable_decomposition::PruningOptions options; // all on, Hom Algorithm A-1/2
options.projected_preimages = false; // opt out of one optimization
stable_decomposition::PruningProfile profile;
auto result = stable_decomposition::pruning(input, epsilon, false, &profile, options);
```

For a preimage, write `F` for the map's generator lift, `S` for the target
submodule generators and `R` for the target presentation. Over F2, form
`D = [F S R]`. A vector `(x,y,z)` belongs to `ker D` exactly when
`F*x = S*y + R*z`; hence the preimage is generated by the **x coordinates** of
the kernel generators. The original implementation constructs all three blocks
of each kernel vector, then discards y and z.

Projected reduction starts the auxiliary matrix at `P = [Id 0 0]`, instead of
the full identity. It applies the same column operations to D and P. If those
operations accumulate to U, the auxiliary matrix is `P*U`; when a column of
`D*U` becomes zero, the corresponding column of `P*U` is its projected kernel
vector. Projection is linear, so these vectors generate exactly the required
preimage; zero projected vectors can be dropped. They need not form a free
kernel basis, and the usual submodule minimization still follows.

This saves tracking, adding and allocating the y/z coefficients, which can grow
during reduction. It does **not** skip reduction of the full D matrix, and it is
not guaranteed to be faster when auxiliary-coordinate work is small.
The implementation reuses the library's existing `bigraded_reduce_pair`.
The [MPfree paper, Section 3 and Algorithms 3/10](https://arxiv.org/html/2010.15623#S3)
describes full-kernel computation with an auxiliary matrix initially equal to
the identity, and the improved grade scheduling in Section 4. This projected
auxiliary construction is a direct adaptation of that invariant; it is not
explicitly presented as a separate algorithm in that paper.

### Composition generators and repeated reductions

Set `N=M(2*epsilon)`, `E=End(M)`, and `H=Hom(M,N)`. The ordinary
`End_d_0` path now computes the actual quotient `Q=H/(eta E)`. It calls
`reduce_matrix_family_modulo(inherited, additional, target_relations)` in
Persistence-Algebra. This normalizes both families modulo the relations
available at each source-generator grade, then performs the usual coefficient
elimination. The two-argument version retains its original coefficient-only
behaviour for ordinary matrix families.

The note distinguishes the default quotient from two subsequent reductions.
`--full-end-space` omits that quotient; it computes actual module maps, still
identifying generator lifts that represent the same map. The optimization
benchmark compares this full-space control with the quotient-based baseline.
Its Quot percentage is `100*(full_space_time/baseline_time - 1)`:
positive values mean that quotienting is faster.

The reductions are:

| Note section | CLI switch | Function |
| --- | --- | --- |
| Appendix A.1 Hom quotient (baseline) | Default; `--full-end-space` disables | `End_d_0` / `reduce_matrix_family_modulo` |
| Composition: Sparse-nilpotent | On; disable with `--no-sparse-nilpotent-generators` | `budgeted_composition_generators` |
| Composition: Complement-orbit | On; disable with `--no-complement-orbit-generators` | `budgeted_composition_generators` |
| Appendix A.3 Adaptive restriction | On; disable with `--no-adaptive-restriction` | `pruning_adaptive_restriction` |

All `End_d_0` overloads take `hom_algorithm` after `info`, using the
library enum `graded_linalg::HomAlgorithm`. The choices follow the library's
function names and [the Hom paper](https://arxiv.org/abs/2604.10381):

| CLI value / enum member | Paper | Library implementation |
| --- | --- | --- |
| `optimised` | A-1/2 | `hom_space_optimised`: restrict Q using local target cokernels |
| `full-restriction` / `full_restriction` | A | `hom_space_full_restriction`: also restrict P using the target's graded kernel |
| `hom-exactness` / `hom_exactness` | B | `hom_space_basis_new(..., true, ...)`: use exactness of Hom |

For example, `End_d_0(eta_d, false, HomAlgorithm::hom_exactness)` selects B.
All three return actual Hom bases and use the same final quotient reduction;
none is claimed to be uniformly faster. Full restriction computes the next
map of the target resolution needed to restrict P, not both entire resolutions.
The input module is `X`, its shifted target is `X_d`, and `d` is the actual
shift: `End_d_0` does not double it. Pruning supplies `d = 2*epsilon`.

Composition selection uses `graded_linalg::budgeted_composition_generators`
and a single prepared quotient-coordinate solver for both sides. It never
constructs the complete pairwise action table.

- **Sparse-nilpotent** samples strictly birth-degree-lowering actions and tries
  small disjoint-support sums. Every certified nilpotent direction can remove an
  original map without retaining the source map.
- **Complement-orbit** tries up to 12 seeds by default, alternating broad sums
  and random subsets of the remaining nonpivot maps. It stops after three
  consecutive low-gain attempts and visits only structurally compatible actions. These
  broader sums can save many more maps, but may cost more in later pruning.
- With both enabled, the second method extends the first method's certificate
  in the original Hom-quotient coordinates. Existing retained sparse candidates
  remain retained; the combined output cannot be larger than the sparse output.

`options.composition_budget.sparse_nilpotent` and
`options.composition_budget.complement_orbit` independently enable the methods;
the same structure exposes their separate work/density budgets. Both default
on. `options.composition_generators=false` disables both. The Hom algorithm and
full-End switches apply to either method, and profiles show both method flags.
The exact complete-action routines remain available for small verification cases.

Search controls are exposed on both the pruning CLI and Hom benchmark:
`--composition-orbit-rounds N` (default 12), `--composition-orbit-seed-mode
sum|random|hybrid` (default hybrid), `--composition-orbit-min-gain N` (default 1),
`--composition-orbit-patience N` (default 3; 0 disables early stopping), and
`--composition-orbit-random-seed N` (default 923407). The default visit allowance
scales with the round limit; `--composition-orbit-pair-visits N` sets an explicit
cap (0 means automatic). Sparse seed counts, nilpotent rounds, product/visit
budgets, and density limits also have `--composition-*` options; see `--help`.
`--all-optimizations` preserves explicitly chosen budgets. Profiles print them.

The optimization experiment additionally compares 3, 6 and 12 hybrid rounds in
`generators.csv` / `generators.tex`. For a wider sweep:

```sh
python3 experiments/benchmark_optimizations.py \
  --composition-orbit-rounds 3 6 12 \
  --composition-orbit-seed-mode sum random hybrid
```

Each sweep configuration enables both composition methods with other optional
improvements off. Search settings and the binary hash are part of the resume key.

The updated local note `explanations/budgeted_composition_generators.pdf` explains
both methods, their shared certificate, complexity, switches, and benchmark
columns. Its editable LaTeX source is alongside it. Earlier general notes and
the manuscript insert remain in `explanations/`.

Total time includes preparation and selection,
which can outweigh the saved iteration work. An empty Q returns `(M, ker eta)`
in both matrix and module pruning, including quick and profiled paths.

Amortized reductions (on by default; `--no-amortize-reductions` disables them)
change K sums only: K stays fixed during a pass, image generators are appended, and the sum is minimized once per batch of 16 images
and at the end. Individual images always skip minimization; this is the
unconditional deferred K image minimization. It removes one reduction per image,
whereas batching reduces how often the entire accumulated sum is minimized.
Batching preserves the exact sum; it trades larger temporary generator matrices
for fewer reductions. It is not a cache of kernel
factorizations across changing preimage problems.

`--adaptive-restriction` computes a smaller family depending on the current pass.
For I and J=eta(I), restrict the chosen maps to `Hom(I,N/J)` and retain a
vector-space basis of their image. For K, do the same in `Hom(K,N/eta(K))`.
It works with a Q basis from either Hom backend or with composition generators. Maps that vanish in these spaces impose no I constraint or add no
new K image modulo eta(K). This computation uses the submodule at the start of
the pass, never its changing intermediate state. Target relations are reduced
with a separate cached scalar basis at each source-generator grade; a later
target relation cannot kill an earlier map value. The profile reports numbers
of candidate and retained restrictions and the time spent selecting them.

`--reuse-target-reductions` changes the witness columns in each I preimage.
At the start of a pass, fix J=eta(I_start). If R is the target relation matrix,
S lifts generators of J, and G_j lifts the current I, then without Reuse each
map f_j computes the projected kernel of `[f_j G_j | S | R]`. The complete
unreduced S and R blocks are supplied again for every map, even though J and
the target module are fixed for the whole pass.

With Reuse, the pass first computes a smaller graded image generating matrix
C' for C=[R | S], preserving the free target's row basis. Each map then uses
`[f_j G_j | C']`. For example, if C has 1000 witness columns and C' has 120,
each preimage uses those 120 target columns. These counts are illustrative;
Reuse pays for one extra target reduction per pass and helps only when the
smaller witness block offsets that cost. It still runs a new kernel reduction
for every changing f_j G_j; it does not cache a kernel factorization. C' is
rebuilt at the next I pass because J may change. K is unaffected. Adaptive
restriction also uses C' when both options are enabled. Projected preimages
reduce auxiliary coordinate tracking, while Reuse reduces the input witness
columns; the two improvements combine.

These new reductions may save map operations while adding preprocessing cost.
Benchmark them individually and in the intended combination; enabling every
switch is not guaranteed to give the fastest run on every input.
The consolidated note supersedes the earlier implementation and composition
notes and includes proofs, examples, and empty-family cases.

### Examples

```bash
# Basic usage with default epsilon
./build/pruning input.scc

# Custom epsilon value
./build/pruning input.scc --epsilon 0.05

# SCC output, paired images, and decomposition comparison
./build/pruning input.scc -e 0.05 --hilbert --aida -o results/pruned.scc

# Only images, without writing another SCC
./build/pruning input.scc -e 0.05 --hilbert --no-output

# Run the original matrix implementation
./build/pruning input.scc -e 0.05 --old

# Compare the two implementations without writing a pruned SCC
./build/pruning input.scc -e 0.05 --compare --no-output
```

## Input Format

The program accepts SCC presentations and projective resolutions supported by
the Persistence-Algebra module reader.

## Output

Output files are automatically named with the pattern: `<input_name>_pru<epsilon>.scc`

Example: `torus3_largestcomp.scc` → `torus3_largestcomp_pru0.0200.scc`

With `-o results/pruned.scc`, optional files use the same directory and stem:

- `pruned_input_hilbert.png` and `pruned_output_hilbert.png`
- `pruned_input_decomposition.sccsum` and `pruned_output_decomposition.sccsum`
- `pruned_decomposition_comparison.txt`
- `pruned_pruning_comparison.txt` with `--compare`

The pruning comparison reports total wall time, I/K iteration counts and average
time per iteration, and call counts, total time and average time for preimages,
intersections, images/compositions, sums, generator reductions, convergence
checks, endomorphism computation, quotient construction, presentation and
minimization. Progress bars remain enabled in both measured runs; phase and total times include their cost.
Input preparation, file I/O, images, AIDA and output checks are outside the timers.
Operation timings include internal library work, so a module intersection's
internal kernel computations belong to its intersection time. The report also
compares generator/relation degree multisets after minimizing both outputs with
the same minimizer. It prints only differing degrees and their multiplicities.
When built with AIDA, it automatically compares the two pruning results' numbers
and types of indecomposables, plus their graded Betti signature multiplicities.
These checks do not write decomposition files. No `--aida` flag is needed;
that flag separately compares the original input with the final pruned module
and saves those decompositions. An exact isomorphism test runs when the maximum
generator-degree multiplicity across both minimized outputs is at most four.
Above that cutoff, the report explicitly marks isomorphism as skipped; matching
invariants alone do not prove isomorphism. Different minimal degree multisets
are rejected immediately regardless of multiplicity. All checks are outside the
pruning timers, and the report is saved even with `--no-output`.

The separate `compare_pruning` regression executable additionally calls the exact
isomorphism test on the minimized outputs, outside the pruning timers, and fails
if they are not isomorphic. Run it with
`ctest --test-dir build -R '^compare_pruning_' --output-on-failure`.
Its CTest cases have a 120-second deadline and always test exact isomorphism.
The ordinary `pruning --compare` uses the multiplicity cutoff because larger
degree blocks can require an exponential search. The cutoff does not bound the
cost of computing Hom on large presentations.

The original matrix `pruning_pair` and `pruning` in `src/pruning.cpp` contain no
profiling code. `src/pruning_profiled.cpp` holds their instrumented copy, used
only for comparisons; algorithm changes must be kept in sync between the two.

The C++ image renderer lives in `Persistence-Algebra/include/grlina/invariants/hilbert_plot.hpp`.
Like `visualisation/visualise_reso.py`, it uses a logarithmic light-blue/blue/black
scale for positive dimensions and white for zero. Both images use the same grid
and colour scale, rounded axis ticks and integer dimension labels. Small sampling
grids are enlarged for readable labels. On macOS, CoreText supplies antialiased
system fonts; other platforms use the built-in bitmap font. No Python or external
plotting library is required. Full projective resolutions are computed explicitly before
evaluating the Hilbert functions, so syzygies are included.

AIDA runs on minimal copies of the input and output. The report compares summand
counts, types, and multiplicities of generator/relation degree signatures at
their actual grades, including pruning's final shift. Matching signatures do
not establish isomorphism; differing decompositions are expected after pruning.

## Default Behavior

An input file is required. Without `--epsilon`, main extracts epsilon from its
presentation. The iteration uses the diagonal shift `(2*epsilon, 2*epsilon)`;
the final result is shifted by `(-epsilon, -epsilon)`. Output filenames record
epsilon, not its doubled value. Epsilon must be finite and nonnegative.

## Version

Current version: **0.3.0** (pre-release)

## Authors

- [Havard Bjerkevik]
- [Jan Jendrysiak](https://github.com/jendjan)
- [Fabian Lenzen](https://gitlab.com/flenzen)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{multiparameter_pruning,
  author = {Havard Bjerkevik and Jan Jendrysiak and Fabian Lenzen},
  title = {Multiparameter Persistence Pruning},
  year = {2025},
  url = {https://github.com/jendjan/Stable-Decomposition},
  license = {GPL-3.0}
}
