
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
- `--fuse-i-preimage` - Compute the I intersection by taking a preimage on I's generator basis.
- `--defer-k-image-minimization` - Minimize each K image when adding it to the sum, rather than also beforehand.
- `--projected-preimages` - Track only required kernel coordinates in two-parameter preimages and intersections.
- `--dense-kernel-pivots` - Use indexed pivot lookup in two-parameter iteration kernels.
- `--throttle-progress` - Print progress when its integer percentage changes.
- `--composition-generators` - Use certified postcomposition generators for I and precomposition generators for K.
- `--amortize-reductions` - Accumulate up to 16 K images before reducing the sum, and skip images of K=0.
- `--hom-classes` - Use actual Hom quotient classes without selecting composition generators.
- `--phase-quotient` - Remove dependent map restrictions separately in each I/K pass.
- `--reuse-target-reductions` - Reduce the fixed I-pass target generators once and reuse them in every preimage.
- `--all-optimizations` - Enable all ten module optimizations.
- `-h, --help` - Show usage and all options.

Options may appear before or after the input. `--epsilon=0.05` also works.

The ten optimization switches default to **off**. Omit a
switch to retain its original behavior, or pass any combination to measure their
interaction. They and `--profile` apply to module pruning and cannot be combined
with `--old` or `--compare`. The legacy matrix algorithms are unchanged.

```bash
# Baseline, then the same module algorithm with one improvement enabled
./build/local-release/pruning input.scc -e 0.05 --no-output --profile
./build/local-release/pruning input.scc -e 0.05 --no-output --profile --fuse-i-preimage

# Paper experiments: run without arguments
python3 experiments/benchmark_runtime.py
python3 experiments/benchmark_optimizations.py
```

Build `build/pruning` separately with CMake in Release mode with compact columns.
Both scripts run that executable and sweep point counts and epsilon percentages on synthetic
and observed real inputs. Runtime uses all optimizations; the optimization
comparison uses an all-off baseline and each of the ten switches enabled alone.
The comparison uses only 80, 160, 320 and 640 points; runtime keeps its wider
20–2048-point sweep.
CSV and LaTeX tables are written automatically, and rerunning continues saved
results. Change the few experiment settings at the top of each Python file.
See [experiments/README.md](experiments/README.md) for data preparation and output paths.

The benchmark uses the module's internal total timer, excluding input preparation
and file I/O. It captures progress output identically for every variant and
shuffles run order with a fixed seed. Thus it includes progress formatting and
log-file writes, but does not measure rendering in an interactive terminal. It
reports timings, not output-equivalence checks. The `pruning_optimizations` CTest
checks all 1024 combinations against the baseline: both submodules of the pruning
pair, iteration counts, and exact isomorphism of the output, including grid and
CSC boundaries. Its executable optionally accepts `FILE EPSILON` to check an
additional input using 32 representative combinations (all eight new-switch
combinations with four earlier settings). Kernel-related switches affect the I/K iteration operations;
the final presentation/minimization retain their existing paths. The composition
switch separately changes construction of the map families.
In other parameter dimensions, the kernel switches retain the existing kernel
implementation.

The C++ API accepts the same options as its last argument:

```cpp
stable_decomposition::PruningOptions options; // all false
options.fuse_i_preimage = true;
options.projected_preimages = true;
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

With `--composition-generators`, set `N=M(2*epsilon)`, `E=End(M)` and
`H=Hom(M,N)`. The ordinary pruning path calls `End_2d_0`, which computes
**generator lifts** modulo transported endomorphism lifts. Wrapping a lift in
`Homomorphism` does not identify lifts that differ by target relations. The
separate `module_hom_space_basis` API already computes actual module-map classes.
The new path uses the existing Hom-exactness algorithm to construct actual H
and E, then forms `Q=H/(eta E)`. Target-relation maps and inherited maps are
reduced in one common coefficient space so that this quotient is well defined.

It computes the action tables `e(2*epsilon)*q` and `q*e` in Q coordinates.
Candidate generators are basis vectors, their total sum and 64 fixed-seed sums.
A greedy search selects candidates whose End orbits increase the span, stopping
only at rank `dim Q`. Basis candidates guarantee termination. This certifies a
generating set, but does not claim the smallest possible one. I uses the left
action (postcomposition); K uses the right action (precomposition). No End maps
are added to the iteration loops. The Hom-class quotient, composition
generators, and pass-specific restriction quotient are explained with examples
and proofs in [composition_generators_pruning.tex](composition_generators_pruning.tex)
and [its PDF](composition_generators_pruning.pdf). The remaining seven
implementation optimizations have their own explanations and proofs in
[pruning_optimizations.tex](pruning_optimizations.tex) and
[its PDF](pruning_optimizations.pdf).

`--hom-classes` stops after constructing a vector-space basis of Q, so its
effect can be measured separately from the End-action search. Composition
selection always includes this class quotient; setting `--hom-classes` as well
does not repeat it. All other switches can be combined independently.

The profile reports the original lift-quotient dimension, actual End and Q
dimensions, and both selected family sizes. Total time includes construction,
action tables and generator selection. Upfront work can outweigh shorter loops;
benchmark the complete algorithm on the intended inputs. The existing
empty-original-family convention is preserved even when Q becomes zero.

`--amortize-reductions` changes K sums only: K stays fixed during a pass, image
generators are appended, and the sum is minimized once per batch of 16 images
and at the end. With `--defer-k-image-minimization`, individual images also skip
minimization. Batching preserves the exact sum; it trades larger temporary
generator matrices for fewer reductions. It is not a cache of kernel
factorizations across changing preimage problems.

`--phase-quotient` computes a smaller family depending on the current pass.
For I and J=eta(I), restrict the chosen maps to `Hom(I,N/J)` and retain a
vector-space basis of their image. For K, do the same in `Hom(K,N/eta(K))`.
It works with either the original lifts, actual Hom classes, or composition
generators. Maps that vanish in these spaces impose no I constraint or add no
new K image modulo eta(K). This computation uses the submodule at the start of
the pass, never its changing intermediate state. Target relations are reduced
with a separate cached scalar basis at each source-generator grade; a later
target relation cannot kill an earlier map value. The profile reports numbers
of candidate and retained restrictions and the time spent selecting them.

`--reuse-target-reductions` forms C=[target relations, eta(I) generators] once
per I pass and computes image generators of C without changing its ambient
row basis. Every preimage then uses that smaller common target matrix. This
shares target preparation, not the complete kernel factorization: each map
still needs its own kernel traversal. The phase quotient also uses the prepared
C when both switches are on. The kernel and projected-kernel switches apply to
these preimages as usual. A proof in the implementation optimization note shows that replacing
C by any matrix with the same graded image preserves every preimage.

These new reductions may save map operations while adding preprocessing cost.
Benchmark them individually and in the intended combination; enabling every
switch is not guaranteed to give the fastest run on every input.
The consolidated note supersedes the earlier individual projection and
composition notes and includes the latest measured examples.

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
