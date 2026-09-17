
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
Use `-DCMAKE_BUILD_TYPE=Debug` for debug symbols; the executable is then named
`pruning_debug`. Build the `tests` and `compare_pruning` targets to run the tests.

The default workspace layout has sibling `Persistence-Algebra`, `AIDA`, and
`Skyscraper-Invariant` checkouts. The build reuses AIDA's bundled CLI11 parser
and Skyscraper-Invariant's `stb_image_write.h`; no Python or GUI library is needed.
Their locations can be overridden with `PERSISTENCE_ALGEBRA_DIR`, `AIDA_DIR`,
`CLI11_INCLUDE_DIR`, and `STB_INCLUDE_DIR`. Boost is required.
`-DPRUNING_WITH_AIDA=OFF` omits the AIDA library integration.

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
- `--image-size <pixels>` - Set the heatmap grid size (default 500; range 128–2048, plus margins for labels).
- `--aida` - Run AIDA directly on both modules and save their decompositions and comparison.
- `--old` - Use the original matrix pruning instead of the module implementation.
- `--quick` - With `--old`, enable the original quick algorithm.
- `--compare` - Run matrix pruning followed by module pruning on independent copies and compare timings and output invariants. AIDA checks run automatically when built in. Both use `quick=false`; the module result is saved. Cannot be combined with `--old`.
- `--no-output` - Skip the pruned SCC; requested images and AIDA results are still saved.
- `-h, --help` - Show usage and all options.

Options may appear before or after the input. `--epsilon=0.05` also works.

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
and saves those decompositions. Matching checks do not prove isomorphism; a
placeholder comment marks where the future isomorphism test belongs. The timing
and comparison report is saved even with `--no-output`.

The original matrix `pruning_pair` and `pruning` in `src/pruning.cpp` contain no
profiling code. `src/pruning_profiled.cpp` holds their instrumented copy, used
only for comparisons; algorithm changes must be kept in sync between the two.

The C++ image renderer lives in `Persistence-Algebra/include/grlina/draw_hf.hpp`.
Like `visualisation/visualise_reso.py`, it uses a logarithmic light-blue/blue/black
scale for positive dimensions and white for zero. Both images use the same grid
and colour scale. Full projective resolutions are computed explicitly before
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

Current version: **0.2**

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
