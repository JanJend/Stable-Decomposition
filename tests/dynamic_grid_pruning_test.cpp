// The shared pruning algorithm must preserve physical grades on coordinate grids.
#define GRLINA_CSC_COMPILE_WARNINGS 0
#define GRLINA_CSC_RUNTIME_WARNINGS 0
#include <doctest/doctest.h>
#include <grlina/mat_csc.hpp>
#include <grlina/mat_grid.hpp>
#include <grlina/mod_epsilon.hpp>
#include <grlina/mod_constructions.hpp>
#include "pruning.hpp"

namespace sd = stable_decomposition;
using namespace graded_linalg;

void compare_grid_pruning() {
    using Matrix = GridGradedSparseMatrix<double, 2, int>;
    using GridModule = Module<Matrix>;
    // Nonuniform coordinates ensure real shifts cannot be implemented as index shifts.
    const sd::Mat free(0, 2, {}, {}, {{0.25, -2}, {1.25, -1}});
    const sd::Mat square(2, 1, {{0}, {0}}, {{0.25, 1.5}, {1.25, -2}}, {{0.25, -2}});
    for (const auto& input : {free, square}) {
        for (double epsilon : {0.0, 0.25, 0.5}) {
            const auto expected = sd::pruning(sd::Module(input), epsilon);
            auto actual = sd::pruning(GridModule(Matrix(input)), epsilon);
            CHECK(actual.presentation().is_graded_matrix());
            const Matrix expected_grid(expected.presentation());
            auto actual_rows = actual.presentation().real_row_degrees();
            auto expected_rows = expected_grid.real_row_degrees();
            auto actual_columns = actual.presentation().real_col_degrees();
            auto expected_columns = expected_grid.real_col_degrees();
            const auto less = Degree_traits<std::array<double, 2>>::lex_lambda();
            std::sort(actual_rows.begin(), actual_rows.end(), less);
            std::sort(expected_rows.begin(), expected_rows.end(), less);
            std::sort(actual_columns.begin(), actual_columns.end(), less);
            std::sort(expected_columns.begin(), expected_columns.end(), less);
            CHECK(actual_rows == expected_rows);
            CHECK(actual_columns == expected_columns);
            for (double x : {-0.25, 0.25, 0.75, 1.25, 1.75, 2.0})
                for (double y : {-2.25, -2.0, -1.5, -1.0, 0.0, 1.5, 2.0})
                    CHECK(actual.dimension_at({x, y}) == expected.dimension_at({x, y}));
        }
    }
}

TEST_CASE("grid pruning agrees with fixed coordinates") {
    compare_grid_pruning();
}

TEST_CASE("module geometry distinguishes independent grids and retains fractional shifts") {
    using Matrix = GridGradedSparseMatrix<double, 2, int>;
    using GridModule = Module<Matrix>;
    Matrix early_matrix(sd::Mat(0, 1, {}, {}, {{0.25, -2}}));
    early_matrix.include_real_degrees({{4.0, 7.0}}); // An unused grid point belongs to the embedding.
    auto early = std::make_shared<const GridModule>(std::move(early_matrix));
    auto late = std::make_shared<const GridModule>(Matrix(sd::Mat(0, 1, {}, {}, {{1.25, -1}})));
    auto sum = direct_sum<Matrix>(early, late);
    CHECK(sum.inclusion_left.check_lifts());
    CHECK(sum.inclusion_right.check_lifts());
    CHECK(sum.projection_left.check_lifts());
    CHECK(sum.projection_right.check_lifts());
    CHECK(sum.module->dimension_at({0.5, -1.5}) == 1);
    CHECK(sum.module->dimension_at({1.25, -1}) == 2);
    auto canonical = Homomorphism<Matrix>::canonical_shift(early, {0.5, 0.5});
    CHECK(canonical.check_lifts());
    CHECK(std::binary_search(canonical.generator_lift().grid(0).begin(), canonical.generator_lift().grid(0).end(), 4.0));
    CHECK(std::binary_search(canonical.target()->presentation().grid(0).begin(), canonical.target()->presentation().grid(0).end(), 3.5));
    CHECK(canonical.generator_lift().real_col_degree(0) == std::array<double, 2>({0.25, -2}));
    CHECK(canonical.generator_lift().real_row_degree(0) == std::array<double, 2>({-0.25, -2.5}));
    auto image = epsilon_image<Matrix>(early, 0.125);
    CHECK(image.generator_map().generator_lift().is_graded_matrix());
    CHECK(image.generator_map().generator_lift().real_col_degree(0) == std::array<double, 2>({0.375, -1.875}));
    CHECK(epsilon_kernel<Matrix>(early, 0.125).is_zero());
    GridModule presented(Matrix(sd::Mat(2, 1, {{0}, {0}}, {{0.25, 1.5}, {1.25, -2}}, {{0.25, -2}})));
    const auto before = presented.hilbert_function({{0.5, -1.5}, {1.25, -1}, {0.5, 1.5}});
    presented.compute_projective_resolution();
    const auto after = presented.hilbert_function({{0.5, -1.5}, {1.25, -1}, {0.5, 1.5}});
    REQUIRE(before.size() == after.size());
    for (std::size_t i = 0; i < before.size(); ++i) CHECK(before[i].dimension == after[i].dimension);
    const auto grid = presented.hilbert_function_on_grid({0.5, 1.25}, {-1.5, 1.5});
    CHECK(grid.values[0][0] == presented.dimension_at({0.5, -1.5}));
    CHECK(grid.values[1][1] == presented.dimension_at({1.25, 1.5}));
}

template <std::size_t Dimension>
void check_dimension_pruning() {
    using Matrix = GridGradedSparseMatrix<double, Dimension, int>;
    using Degree = typename Matrix::real_degree_type;
    Degree early{}, late{}, expected{}, below{};
    early.fill(0.25); late.fill(1.25); expected.fill(1.75); below.fill(1.5);
    Matrix matrix(0, 2, {}, {early, late});
    auto result = sd::pruning(Module<Matrix>(std::move(matrix)), 0.5);
    CHECK(result.presentation().is_graded_matrix());
    CHECK(result.number_of_generators() == 2);
    CHECK(result.number_of_relations() == 0);
    for (const auto& degree : result.presentation().real_row_degrees())
        CHECK(degree == expected);
    CHECK(result.dimension_at(below) == 0);
    CHECK(result.dimension_at(expected) == 2);
}

TEST_CASE("pruning uses the diagonal in one and three parameters") {
    check_dimension_pruning<1>();
    check_dimension_pruning<3>();
}

TEST_CASE("grid pruning retains the scalar and index template choices") {
    using Matrix = GridGradedSparseMatrix<float, 2, long>;
    Matrix matrix(0, 2);
    matrix.set_real_degrees({}, {{0.25f, -2.0f}, {1.25f, -1.0f}});
    auto result = sd::pruning(Module<Matrix>(std::move(matrix)), 0.5);
    CHECK(result.presentation().is_graded_matrix());
    CHECK(result.number_of_generators() == 2);
    CHECK(result.dimension_at({1.5f, -0.5f}) == 0);
    CHECK(result.dimension_at({1.75f, -0.5f}) == 2);
    const auto grid = result.hilbert_function_on_grid({1.5, 1.75}, {-0.5, 0.0});
    CHECK(grid.values[0][0] == 0);
    CHECK(grid.values[1][0] == 2);
}
