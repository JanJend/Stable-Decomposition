#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <string>
#include "../include/pruning.hpp"

TEST_CASE("test1.scc") {
    using namespace stable_decomposition;
    std::string input = "test1.scc";
    float delta = .01;

    // The module overload must be an owning adapter for the established
    // matrix algorithm; compare both paths directly rather than relying on an
    // older generated fixture whose degree-shift convention became stale.
    Mat matrix_input(input);
    Mat matrix_result = pruning(matrix_input, delta);
    PModule module_result = pruning(PModule(input), delta);
    CHECK(module_result.presentation().row_degrees == matrix_result.row_degrees);
    CHECK(module_result.presentation().col_degrees == matrix_result.col_degrees);
    CHECK(module_result.presentation().data == matrix_result.data);

    auto parent = std::make_shared<PModule>(input);
    auto [I, K] = pruning_pair(parent, delta);
    CHECK(I.parent() == parent);
    CHECK(K.parent() == parent);
    CHECK(I.generators().row_degrees == parent->presentation().row_degrees);
    CHECK(K.generators().row_degrees == parent->presentation().row_degrees);

    //TODO F: One cannot really check if P and P_expected present isomorphic modules, because that's a hard problem.
    //        Testing if two matrices span /identical/ submodules or subquotients of a given module should be easier;
    //        for example, one can do:
    // 
    // auto [I, K] = pruning_pair(M, delta);
    // Mat I_expected(...), K_expected(...);
    // CHECK(present_same_submodule(M, I, I_expected));
    // CHECK(present_same_submodule(M, K, K_expected));

}
