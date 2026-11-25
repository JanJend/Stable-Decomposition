# Multiparameter Persistence Pruning

A C++ program for pruning two-parameter persistence modules. It reads an SCC
presentation or projective resolution and writes a minimal SCC presentation of
the pruned module. It can also plot Hilbert functions and compare decompositions
using AIDA.

## Build

You need CMake 3.28 or newer, a C++17 compiler, Make, and Boost (timer, graph,
and chrono). Keep the dependency checkouts next to this repository:

```text
workspace/
├── Stable-Decomposition/
├── Persistence-Algebra/
└── AIDA/
```

AIDA supplies the CLI11 argument parser. Persistence-Algebra supplies the algebra
and PNG writer. You do not need to build either dependency separately.

From the `Stable-Decomposition` directory:

```sh
mkdir -p build
cd build
cmake ..
make
cd ..
```

- `mkdir -p build` creates a directory for compiled files; it also works if the directory already exists.
- `cmake ..` reads `CMakeLists.txt` from the parent directory, finds dependencies, and generates the Makefiles in `build/`.
- `make` compiles the program and enabled tests.
- `cd ..` returns to the repository root, where the examples below run.

The executable is `build/pruning`. The default build already uses Release mode,
compact matrix columns, and AIDA support. No extra CMake options are needed.
After changing C++ code, run `make` inside `build/` again.

## Run

For example, using an included input:

```sh
./build/pruning test_presentations/test1.scc
```

Choose epsilon as a percentage of the presentation's bounding box, or supply
its coordinates directly:

```sh
./build/pruning input.scc --epsilon-percent 1
./build/pruning input.scc --epsilon-vector 0.05 0.002 -o results/pruned.scc
```

`--epsilon-percent P` computes `(P/100 * Lx, P/100 * Ly)`, where `Lx` and `Ly`
are the bounding-box side lengths across all generator and relation degrees.
The default is 1%. A zero-length side gives a zero component; an empty
presentation gives `(0,0)`.

`--epsilon-vector X Y` supplies absolute amounts in the input's coordinate
units. The existing `--epsilon E` (or `--delta E`) is a shortcut for `(E,E)`.
These three input modes are mutually exclusive. Every component and percentage
must be finite and nonnegative, and the doubled shift must be representable.

Without `-o`, the result is saved beside the input as
`<input_name>_pru<X>_<Y>.scc`, with components written to four decimal places.
Equal components retain the shorter `_pru<E>.scc` name. Missing output
directories are created automatically.

The algorithm uses the componentwise shift `2 * epsilon` and translates
the final result by `-epsilon`. For a box of size 100 by 10, 1% gives epsilon
`(1,0.1)` and the internal shift `(2,0.2)`.

Pruning automatically detects the available logical CPU cores at startup and
uses them as its worker budget. If detection fails, it uses one worker. This
works with the default build; no extra CMake settings or OpenMP are needed.
Use `--threads 1` for a serial run, or `--threads N` to choose another limit.
With `--profile`, I/K phase times can overlap and should not be added together.

### Optional outputs and controls

Add these only when you need them:

| Option | Use it to… |
| --- | --- |
| `--hilbert` | Save input and output Hilbert PNGs with shared axes and colour scale. |
| `--image-size 1000` | Change the Hilbert sampling grid size; default 500, allowed range 128–2048. |
| `--aida` | Decompose the input and output and save a comparison report. |
| `--no-output` | Skip the pruned SCC while keeping requested images and reports. |
| `--profile` | Print operation timings and optimization settings. |
| `--threads 4` | Override automatic core detection with a limit of four computation workers. |
| `--memory-limit-mb 4096` | Stop if sampled memory use exceeds 4096 MiB; default 0 disables the limit. |
| `--old` | Run the original matrix implementation; add `--quick` for its quick mode. |
| `--compare` | Run both implementations and save timings and output checks; save the module implementation's result. |
| `--help` | Show all options. |

For SCC output, plots, and AIDA decompositions together:

```sh
./build/pruning input.scc -e 0.05 --hilbert --aida -o results/pruned.scc
```

This writes `results/pruned.scc` and, in the same directory:

- `pruned_input_hilbert.png` and `pruned_output_hilbert.png`
- `pruned_input_decomposition.sccsum` and `pruned_output_decomposition.sccsum`
- `pruned_decomposition_comparison.txt`

AIDA compares summand counts, types, and graded Betti signatures. Matching
signatures do not prove isomorphism. `--compare` writes a separate
`*_pruning_comparison.txt`; it reports when an exact isomorphism check is skipped.

### Algorithm experiments

The module implementation enables its optional optimizations by default, except fast I, and
uses Hom Algorithm A-1/2 (`optimised`). It computes the quotient with the
single-solve `End_d_0` routine, while retaining a separate End₀ basis for the
left and right composition actions. To compare algorithms or measure an
optimization, use `--hom-algorithm hom-exactness` (B),
`--hom-algorithm full-restriction` (A), or the `--no-...` switches listed by
`--help`. Composition search budgets are also listed there.

The tuned **Nilpotent reduction** and **Complement-orbit** searches are both
enabled by default. Nilpotent reduction first certifies individual quotient
basis maps in the image of the nilpotent ideal using elementary factorizations.
It tries known valid factors before checking new ones against the relations.
Successful certificates insert the maps' known quotient coordinates directly.
Sampled lowering products use the remaining allowance; small ordinary trials
are skipped when direct certificates fill that allowance. Complement-orbit
then tries broad sums and random subsets with its own limits.

The full quotient-coordinate solver is prepared only when a product or a
new-factor membership check needs it, and is shared by both sides. With the
ordinary defaults, direct certificates and sampled products share a 4,096-unit
allowance per side; certificates plus the small-trial products share an
8,192-unit total. The certificate search has a separate 8,192-visit cap.
Every returned family still generates the full quotient over End₀.
See the [exact production algorithm](explanations/composition_generators_implementation.tex)
and [detailed explanation](explanations/budgeted_composition_generators.pdf)
for the checks, candidate order, budgets, and correctness proof.

Use `--no-single-solve-quotient` to compare with the separate End₀ and shifted-Hom
computation; the optimization experiment reports this comparison as **Direct Q**.

The I loop uses the **standard** rule by default: fix
`J = shift(I)` at the start of each pass. Use `--update-i-target` for the
**fast** rule, refreshing `J` from the current I after its embedded generators
change, or
`--no-update-i-target` for standard. `--all-optimizations` also enables fast I.
The optimization experiment compares it with the standard baseline in the
**Fast I** column; K uses the same update rule in both modes.

Pruning keeps map families in sparse coefficient form and constructs temporary
lifts on demand. Adaptive restriction selection processes bounded batches and
retains sparse independence pivots, avoiding a full family of restricted matrices.

`--profile`, `--threads`, and algorithm tuning options apply to the module
implementation and cannot be combined with `--old` or `--compare`.
See [the experiment guide](experiments/README.md) for benchmark setup and runs.

## Tests

Tests are built by default. After `make`, run them from the repository root:

```sh
cd build
ctest --output-on-failure
cd ..
```

`ctest` runs the tests; `--output-on-failure` prints details for failures.
Installing doctest enables additional tests. The standalone regression tests
work without it. See [functional tests](tests/functional/README.md) and the
[hand-checked examples](test_presentations/tests_with_solutions/README.md).

## Optional build settings

Use these only for a specific build requirement. Inside `build/`, add the
setting to `cmake ..`, then run `make`. For example, to build with debug symbols:

```sh
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
cd ..
```

This creates `build/pruning_debug`; use that path when running the debug build.

`-DNAME=value` sets a CMake configuration value. CMake remembers these settings
in `build/CMakeCache.txt`, so a later `cmake ..` keeps them. Change a value
explicitly to switch back.

| Setting | When to use it |
| --- | --- |
| `-DCMAKE_BUILD_TYPE=Debug` | Debug with symbols and checks. Set `Release` to return to the default. |
| `-DBUILD_TESTING=OFF` | Skip building tests. Set `ON` to restore them. |
| `-DPRUNING_WITH_AIDA=OFF` | Omit decomposition support. CLI11 is still required; set `ON` to restore AIDA. |
| `-DGRLINA_USE_COMPACT_COLUMNS=OFF` | Compare ordinary vector columns with the default compact storage. Set `ON` to restore compact columns. |
| `-DPERSISTENCE_ALGEBRA_DIR=/path/to/Persistence-Algebra` | Use an algebra checkout outside the sibling layout. |
| `-DAIDA_DIR=/path/to/AIDA` | Use an AIDA checkout outside the sibling layout. |
| `-DCLI11_INCLUDE_DIR=/path/to/include` | Supply the directory containing `CLI11.hpp` separately. |
| `-DBoost_ROOT=/path/to/boost` | Select a particular Boost installation if automatic discovery chooses the wrong one. |

For faster compilation, use `make -j4` to run up to four compile jobs. To build
only the command-line program, use `make pruning`. These commands run inside
`build/`. The `--threads` runtime option controls pruning workers separately.

On Apple Silicon, Boost must match the compiler's ARM architecture. If CMake
has cached the wrong Boost installation, add `-U '*[Bb]oost*'` alongside
`-DBoost_ROOT=...` to clear those cached paths. A `Could NOT find OpenMP` message
is harmless: OpenMP is optional.

## Further reading

- [Experiments](experiments/README.md): data generation, runtime measurements, optimization comparisons, and Hom-size benchmarks.
- [Functional tests](tests/functional/README.md): implementation comparisons and regression checks.
- [Algebra extraction](docs/algebra-extraction.md): the split between this project and Persistence-Algebra.

`main.cpp` shows the pruning workflow. Implementations live in `src/`, headers
in `include/`, and example inputs in `test_presentations/`.

## Authors and license

Havard Bjerkevik, [Jan Jendrysiak](https://github.com/jendjan), and
[Fabian Lenzen](https://gitlab.com/flenzen).

Licensed under the [GNU General Public License v3.0 or later](LICENSE).

## Citation

```bibtex
@software{multiparameter_pruning,
  author = {Havard Bjerkevik and Jan Jendrysiak and Fabian Lenzen},
  title = {Multiparameter Persistence Pruning},
  year = {2025},
  url = {https://github.com/jendjan/Stable-Decomposition},
  license = {GPL-3.0}
}
```
