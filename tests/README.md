# Comparing the two pruning implementations

From the Stable-Decomposition directory:

```sh
cmake -S . -B build/vscode-debug -DBUILD_TESTING=ON
cmake --build build/vscode-debug --target compare_pruning -j 2
./build/vscode-debug/compare_pruning
./build/vscode-debug/compare_pruning --epsilon 0.25 toy_example_2.scc
./build/vscode-debug/compare_pruning --epsilon 0.01 tests/test1.scc tests/test2.scc
ctest --test-dir build/vscode-debug -R '^compare_pruning_' --output-on-failure
```

With no filenames, the program uses `toy_example_1.scc`, `toy_example_2.scc`,
and `no_columns_test.scc` from the configured Persistence-Algebra checkout's
`test_presentations` directory. Bare filenames are also resolved there if they
do not exist in the current directory. Explicit paths can point to any input
accepted by the module reader. The default epsilon is **0.5**; the default is
the ordinary matrix algorithm (`quick=false`). `--quick` passes `true` to both
overloads, although the module implementation currently ignores this flag.

Both overloads receive independent copies of the same sorted presentation.
The program validates their outputs, reports runtimes and exact matrix equality,
then compares:

- Generator and relation counts and degree multisets after minimizing copies
  with the same minimizer (graded Betti degrees in homological degrees 0 and 1).
- Fibre dimensions of the original outputs on the Cartesian product of their
  combined critical coordinates, including points before and after the extrema.
  For finite presentation degrees, the full critical grid covers the regions
  where the Hilbert function can change. Axes exceeding 64 coordinates are
  sampled evenly and the report explicitly labels the grid as sampled.

Degree comparisons use exact doubles, so numerical discrepancies are visible.
Up to eight differing degrees and eight fibre mismatches are printed per input.
Different matrix entries alone do not fail the comparison: bases can differ.
Matching invariants **do not prove module isomorphism**, even on the full grid.
No presentation files are modified or written.

Exit codes: **0** for matching checked invariants, **1** for a mismatch, **2** for
invalid arguments or a caught runtime error. Multiple inputs are all checked
unless an underlying library assertion or abort terminates the process.

The `pruning_cli_outputs` CTest exercises SCC round trips, PNG output, CLI help
and option errors, and (when enabled) direct AIDA decomposition and comparison.
It covers two free summands and the zero module, writing fixtures under the
build directory. Build the `pruning` target before running this test.
