# Algebra extraction and source layout

## Layout

`main.cpp` remains the executable entry point at the repository root. This is a
conventional small-project layout, not a C++ requirement. All reusable compiled
implementation now lives in `src/`; public declarations remain in `include/`.

- `src/pruning.cpp`: only pruning-pair iterations and pruning entry points.
- `src/delta.cpp` / `include/delta.hpp`: the parameter-selection heuristic.
- `src/utils.cpp` / `include/utils.hpp`: command-line parsing and output filenames.
- `src/progress.cpp` / `include/progress.hpp`: iteration-progress display.
- `src/algebra_compat.cpp` / `include/algebra_compat.hpp`: thin wrappers retaining
  the historical Stable-Decomposition algebra names.
- `include/types.hpp`: the shared matrix/module aliases.

CMake compiles these sources once into `stable_pruning`. Both `pruning` and
`tests` link that library. The VS Code launch still builds the `pruning` target;
it follows the source moves automatically. Set algorithm breakpoints in
`src/pruning.cpp`, not the removed root-level file. No nested library checkout
was removed or updated.

## General operations now owned by Persistence-Algebra

| Historical helper | Library implementation / preferred object API |
| --- | --- |
| `matrix_reduction` | `reduce_matrix_family_modulo` in `matrix_family.hpp` |
| `homSpace` | `homomorphism_lift_basis`; `module_hom_space_basis` for Hom classes |
| `End_2d_0` | `shifted_endomorphism_lift_complement` |
| `zero_submodule`, `all_submodule` | matrix adapters; `Submodule::zero`, `Submodule::whole` |
| `submodule_sum` | free-target matrix adapter; `Submodule::sum` for a presented parent |
| `reduce_submodule` | adapter to `Submodule::minimize_generators` |
| `shifting_morphism` | `canonical_shift_lift`; `Homomorphism::canonical_shift` |
| shifting an existing map | `Homomorphism::shifted`, including all stored lifts |
| `image_contained_in_image` | exact degree-admissible linear system for a free target |
| `present_same_submodule` | `Submodule::equals`, including parent relations |
| containment in a presented parent | `Submodule::contains`, `is_contained_in` |
| `image(f,A,B,U)` (previously only declared) | `Homomorphism::image(submodule)` |
| `timed_with_progress` | shared exception-safe implementation in `grlina/progress.hpp` |

Presentation-level adapters live in `grlina/presentation_operations.hpp`, also
available through `grlina/modules.hpp`. Object operations live in `submodule.hpp`
and `homomorphism.hpp`. All matrix templates require the graded CRTP contract.
Old Stable-Decomposition function signatures remain available through
`pruning.hpp` and link via `stable_pruning`; their mathematical implementations
are no longer duplicated here.

## Mathematical distinctions and behavioral changes

Containment in a presented parent P solves `[P B] X = A`, with each solution
column restricted to degrees <= the corresponding degree of A. Equality tests
both inclusions. Free-target inclusion instead solves `B X = A`. These checks
require neither sorting nor a graded kernel. Parent-aware operations require
the same parent object to fix the ambient basis.

The old containment test compared kernel generator degree multisets. It is
replaced by the direct membership test above. The free-target sum adapter also
preserves the original ambient row basis; it no longer sorts that basis silently.

Matrix-family reduction treats matrices as coefficient vectors, not graded
columns. It preserves span(A), reduces B modulo A, and returns independent
representatives. Pivot positions are pairs of indices, so rectangular matrices
do not use the old mismatched flatten/unflatten convention. Reducers reference
stable input storage; surviving B matrices are moved only after reduction.

The shifted-lift complement retains the existing pruning convention: it takes
a quotient of spaces of generator lifts, not Hom classes modulo target
relations. Unshifted lifts are regraded into the shifted target before reduction.
The old `End_2d_0` name is retained only for compatibility; its actual shift is
`(delta, delta)`. This extraction does not change that convention or claim to
prove the pruning algorithm's mathematical correctness.

Canonical shift maps subtract the shift from target generator degrees. They
include identity lifts at every stored projective group, reject inadmissible
shifts, and retain the module's resolution. Shifting an existing homomorphism
translates both modules and all supplied lifts without modifying the originals.

Quick pruning now returns `(whole, zero)` immediately when the additional lift
family is empty; previously its loop could never reach its stopping condition.
The timer joins its worker on exceptions and wakes immediately on completion;
it no longer delays every fast operation by a polling interval.

## Delta behavior deliberately preserved

The extracted heuristic is unchanged: use the explicit delta if supplied;
otherwise take 1% of the maximum coordinate extent over row and column degrees,
falling back to 0.01 for empty/degenerate inputs. **main.cpp still explicitly
sets delta to 1 and leaves the get_delta call commented out.** This refactor
does not override that existing experimental choice; using the heuristic or
`--delta` in the executable requires restoring that call.

## Tests

`Persistence-Algebra/tests/presentation_operations_test.cpp` has handmade cases
for incomparable generator degrees, equality only modulo relations, unsorted
ambient bases, no-kernel R4 containment, rectangular matrix families, shifted
lifts, retained resolutions and restricted images. The existing module tests
continue to cover the underlying kernel/minimization machinery.

Stable-Decomposition's `algebra_adapters_test.cpp` exercises every old wrapper,
zero-scale quick pruning, and timer success/reference/void/move-only/exception
paths. `delta_test.cpp` checks empty, degenerate, row-only and mixed-degree
heuristics and explicit overrides. Existing pruning integration tests remain.
