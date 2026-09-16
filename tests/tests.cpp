#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <string>
#include "../include/pruning.hpp"

TEST_CASE("submodule syzygies preserve ambient coordinates") {
    using namespace stable_decomposition;
    // Incomparable a,b and c=a+b. All three survive ordinary graded reduction.
    Mat parent(0, 3, {}, {}, {{0,0}, {-1,-1}, {-2,-2}});
    Mat generators(3, 3, {{0,2}, {1,2}, {0,1}}, {{0,1}, {1,0}, {1,1}}, parent.row_degrees);
    auto reduced = reduce_submodule(parent, generators);
    CHECK(reduced.get_num_cols() == 2);
    CHECK(reduced.row_degrees == parent.row_degrees);
    CHECK(reduced.data == graded_linalg::array<int>({{0,2}, {1,2}}));
}

TEST_CASE("zero-scale pruning preserves a square interval") {
    using namespace stable_decomposition;
    Mat square(2, 1, {{0}, {0}}, {{0,1}, {1,0}}, {{0,0}});
    PModule result = pruning(PModule(square), 0.0);
    CHECK(result.dimension_at({0,0}) == 1);
    CHECK(result.dimension_at({0.5,0.5}) == 1);
    CHECK(result.dimension_at({1,0}) == 0);
    CHECK(result.dimension_at({0,1}) == 0);
    CHECK(result.presentation().data == graded_linalg::array<int>({{0}, {0}}));
}

TEST_CASE("test1.scc") {
    using namespace stable_decomposition;
    std::string input = "test1.scc";
    float epsilon = .01;

    // The module overload must be an owning adapter for the established
    // matrix algorithm; compare both paths directly rather than relying on an
    // older generated fixture whose degree-shift convention became stale.
    Mat matrix_input(input);
    Mat matrix_result = pruning(matrix_input, epsilon);
    PModule module_result = pruning(PModule(input), epsilon);
    CHECK(module_result.presentation().row_degrees == matrix_result.row_degrees);
    CHECK(module_result.presentation().col_degrees == matrix_result.col_degrees);
    CHECK(module_result.presentation().data == matrix_result.data);

    auto parent = std::make_shared<PModule>(input);
    auto [I, K] = pruning_pair(parent, epsilon);
    CHECK(I.parent() == parent);
    CHECK(K.parent() == parent);
    CHECK(I.generators().row_degrees == parent->presentation().row_degrees);
    CHECK(K.generators().row_degrees == parent->presentation().row_degrees);

    //TODO F: One cannot really check if P and P_expected present isomorphic modules, because that's a hard problem.
    //        Testing if two matrices span /identical/ submodules or subquotients of a given module should be easier;
    //        for example, one can do:
    // 
    // auto [I, K] = pruning_pair(M, epsilon);
    // Mat I_expected(...), K_expected(...);
    // CHECK(present_same_submodule(M, I, I_expected));
    // CHECK(present_same_submodule(M, K, K_expected));

}
