This folder contains presentations and prunings (up to isomorphism) of the examples with a given pruning parameter.
Note that different presentations may give rise to the same module up to iso.

The solutions with name [filename]_pru[pruning parameter].scc where the pruning parameter does not have trailing zeros have been checked by hand. (Which still leaves the possibility that they're wrong.)

test4.scc is a version of the example from the original pruning paper.

The historical `_pru1` suffix records the iteration shift, which is
`2 * epsilon` in the current API. Run these examples with `epsilon = 0.5`;
their final translation is `(-0.5, -0.5)`.

`test1_pru1.scc` is the zero module. Its formerly empty file now has the SCC
preamble and zero dimensions so it can be read like the other solutions.

The regression in `tests/pruning_solutions_test.cpp` computes all four prunings
and compares them with these solutions using the exact isomorphism test.
