#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <string>
#include "../include/pruning.hpp"

TEST_CASE("test1.scc") {
    using namespace stable_decomposition;
    std::string input = "test1.scc", result = "test1_pru0.0100.scc";
    float delta = .01;
    Mat M(input);
    auto P = pruning(M, delta);
    Mat P_expected(result);
    CHECK(P.row_degrees == P_expected.row_degrees && P.col_degrees == P_expected.col_degrees);

    //TODO F: One cannot really check if P and P_expected present isomorphic modules, because that's a hard problem.
    //        Testing if two matrices span /identical/ submodules or subquotients of a given module should be easier;
    //        for example, one can do:
    // 
    // auto [I, K] = pruning_pair(M, delta);
    // Mat I_expected(...), K_expected(...);
    // CHECK(present_same_submodule(M, I, I_expected));
    // CHECK(present_same_submodule(M, K, K_expected));

}
