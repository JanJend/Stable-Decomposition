// Independent switches must preserve the pruning pair, output and iteration counts.
#define GRLINA_CSC_COMPILE_WARNINGS 0
#define GRLINA_CSC_RUNTIME_WARNINGS 0
#include "pruning.hpp"
#include <grlina/isomorphism_test.hpp>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>

namespace sd = stable_decomposition;
namespace gl = graded_linalg;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

sd::PruningOptions options_for(unsigned mask) {
    return {bool(mask & 1), bool(mask & 2), bool(mask & 4), bool(mask & 8), bool(mask & 16),
            bool(mask & 32), bool(mask & 64), bool(mask & 128), bool(mask & 256), bool(mask & 512)};
}

template<class Matrix>
sd::Mat coordinate_matrix(const Matrix& matrix) {
    gl::array<int> columns;
    for (int j = 0; j < matrix.get_num_cols(); ++j) columns.push_back(matrix.get_col(j));
    return sd::Mat(matrix.get_num_cols(), matrix.get_num_rows(), columns,
                   gl::geometric_col_degrees(matrix), gl::geometric_row_degrees(matrix));
}

template<class Matrix>
void check_pruning(const Matrix& matrix, double epsilon, bool exhaustive = true) {
    using Module = gl::Module<Matrix>;
    auto parent = std::make_shared<const Module>(matrix);
    sd::PruningProfile baseline_profile;
    const auto baseline_pair = sd::pruning_pair(*parent, epsilon, false, &baseline_profile);
    const auto baseline = sd::pruning(Module(matrix), epsilon);
    std::vector<unsigned> masks;
    if (exhaustive) {
        for (unsigned mask = 0; mask < 1024; ++mask) masks.push_back(mask);
    } else {
        // Larger inputs cover all eight new combinations with four prior
        // settings; small and random fixtures cover all 1024 combinations.
        for (unsigned added = 0; added < 8; ++added)
            for (unsigned prior : {0u, 13u, 32u, 127u}) masks.push_back(prior | (added << 7));
    }
    for (unsigned mask : masks) {
        const auto options = options_for(mask);
        sd::PruningProfile profile;
        const auto pair = sd::pruning_pair(*parent, epsilon, false, &profile, options);
        require(pair.first.equals(baseline_pair.first), "I changed under an optimization");
        require(pair.second.equals(baseline_pair.second), "K changed under an optimization");
        require(profile.i_iterations == baseline_profile.i_iterations &&
                profile.k_iterations == baseline_profile.k_iterations, "Iteration counts changed");
        require(profile.i_phase_kept <= profile.i_phase_candidates &&
                profile.k_phase_kept <= profile.k_phase_candidates, "Phase quotient increased its family");
        const auto result = sd::pruning(Module(matrix), epsilon, false, nullptr, options);
        require(result.presentation().is_graded_matrix(), "Ungraded optimized result");
        require(gl::is_isomorphic(coordinate_matrix(baseline.presentation()),
                                   coordinate_matrix(result.presentation())), "Pruning output changed");
    }
}

void check_projected_preimages() {
    std::mt19937 random(42);
    for (int example = 0; example < 40; ++example) {
        const int rows = example % 6;
        gl::vec<gl::r2degree> row_degrees;
        for (int r = 0; r < rows; ++r)
            row_degrees.push_back({double(random() % 4), double(random() % 4)});
        auto make_matrix = [&](int count) {
            gl::vec<gl::r2degree> degrees;
            gl::array<int> columns(count);
            for (int c = 0; c < count; ++c) {
                degrees.push_back({double(random() % 5), double(random() % 5)});
                for (int r = 0; r < rows; ++r)
                    if (gl::Degree_traits<gl::r2degree>::smaller_equal(row_degrees[r], degrees.back()) && random() % 2)
                        columns[c].push_back(r);
            }
            return sd::Mat(count, rows, columns, degrees, row_degrees);
        };
        auto map = make_matrix(example % 7);
        auto presentation = make_matrix(example % 5);
        auto submodule = make_matrix(example % 4);
        auto domain = std::make_shared<const sd::Module>(
            sd::Mat(0, map.get_num_cols(), {}, {}, map.col_degrees));
        sd::Submodule expected(domain, map.inverse_image_copy(presentation, submodule));
        for (bool projected : {false, true}) for (bool dense : {false, true}) {
            sd::Submodule actual(domain, gl::inverse_image_generators(map, presentation, submodule, projected, dense));
            require(actual.equals(expected), "Projected/dense preimage differs from full kernel");
            require(actual.generator_map().generator_lift().row_degrees == map.col_degrees,
                    "Preimage changed its domain basis");
            auto combined = presentation;
            combined.append_matrix(submodule);
            auto reduced = sd::pruning_reduce_target(std::move(combined), dense);
            require(reduced.row_degrees == map.row_degrees, "Target preparation changed the ambient basis");
            sd::Submodule reused(domain, sd::pruning_preimage_with_target(map, reduced, projected, dense));
            require(reused.equals(expected), "Reused target changed the preimage");
        }
    }
}

void check_phase_quotient_grades() {
    auto parent = std::make_shared<const sd::Module>(sd::Mat(0, 1, {}, {}, {{0, 0}}));
    sd::Mat identity(1, 1, {{0}}, {{0, 0}}, {{0, 0}});
    std::vector<sd::Hom> maps{sd::Hom(parent, parent, identity), sd::Hom(parent, parent, identity)};
    sd::Mat target(1, 1, {{0}}, {{1, 1}}, {{0, 0}});
    require(sd::pruning_phase_map_basis(maps, identity, target) == std::vector<std::size_t>{0},
            "A later target generator incorrectly killed an earlier map restriction");
    sd::Mat later(1, 1, {{0}}, {{2, 2}}, {{0, 0}});
    require(sd::pruning_phase_map_basis(maps, later, target).empty(),
            "A target generator did not kill a restriction at a later grade");
}

template<std::size_t Dimension>
void check_other_dimensions() {
    using Matrix = gl::GridGradedSparseMatrix<double, Dimension, int>;
    using Module = gl::Module<Matrix>;
    typename Matrix::GridD early{}, late{};
    early.fill(0.25); late.fill(1.25);
    auto parent = std::make_shared<const Module>(Matrix(0, 2, {}, {early, late}));
    sd::PruningProfile expected_profile;
    const auto expected = sd::pruning_pair(*parent, 0.5, false, &expected_profile);
    for (unsigned added = 0; added < 8; ++added) for (unsigned prior : {0u, 127u}) {
        sd::PruningProfile profile;
        const auto result = sd::pruning_pair(*parent, 0.5, false, &profile, options_for(prior | (added << 7)));
        require(result.first.equals(expected.first) && result.second.equals(expected.second),
                "New switches changed pruning in another coordinate dimension");
        require(profile.i_iterations == expected_profile.i_iterations && profile.k_iterations == expected_profile.k_iterations,
                "New switches changed iteration counts in another coordinate dimension");
    }
}

void check_composition_generators() {
    // Use a small artificial conversion limit to exercise the uncapped path.
    sd::Mat free(0, 2, {}, {}, {{0, 0}, {0, 0}});
    free.compute_rows_forward();
    require(gl::hom_space_basis_new(free, free, true, false, true, 1).empty(),
            "Explicit Hom-space conversion limit was ignored");
    require(gl::hom_space_basis_new(free, free, true, false, true, 0).size() == 4,
            "Uncapped Hom-space basis lost endomorphisms");
    // Exercise quotient coordinates beyond a machine word: 10 incomparable
    // free generators give dim Hom=100, dim End=10, hence dim Q=90.
    gl::vec<gl::r2degree> degrees;
    for (int j = 0; j < 10; ++j) degrees.push_back({double(j), double(9-j)});
    auto parent = std::make_shared<const sd::Module>(sd::Mat(0, 10, {}, {}, degrees));
    const auto eta = sd::Hom::canonical_shift(parent, {10, 10});
    const auto families = sd::pruning_composition_generators(eta, gl::r2degree{10, 10});
    require(families.original_size == 90 && families.quotient_dimension == 90,
            "Incorrect large Hom quotient dimension");
    require(!families.i_maps.empty() && families.i_maps.size() < 90 &&
            !families.k_maps.empty() && families.k_maps.size() < 90,
            "Composition search did not reduce the nontrivial map families");

    std::mt19937 random(20261103);
    for (int example = 0; example < 12; ++example) {
        gl::vec<gl::r2degree> rows, columns;
        gl::array<int> relations(3);
        for (int r = 0; r < 3; ++r)
            rows.push_back({double(random() % 3), double(random() % 3)});
        for (int c = 0; c < 3; ++c) {
            columns.push_back({double(random() % 4), double(random() % 4)});
            for (int r = 0; r < 3; ++r)
                if (gl::Degree_traits<gl::r2degree>::smaller_equal(rows[r], columns[c]) && random() % 2)
                    relations[c].push_back(r);
        }
        check_pruning(sd::Mat(3, 3, relations, columns, rows), double(example % 4) / 4);
    }
}

int main(int argc, char** argv) {
    try {
        check_projected_preimages();
        check_phase_quotient_grades();
        check_other_dimensions<1>();
        check_other_dimensions<3>();
        check_composition_generators();
        const std::filesystem::path fixtures(PRUNING_TEST_FIXTURE_DIR);
        for (int i = 1; i <= 4; ++i) {
            const auto path = fixtures / "tests_with_solutions" / ("test" + std::to_string(i) + ".scc");
            check_pruning(sd::Module(path.string()).presentation(), 0.5);
        }
        const sd::Module nonzero_k(std::string(PRUNING_TEST_PRESENTATIONS_DIR) + "/small_pres_with_nonzero_k_min_pres.firep");
        check_pruning(nonzero_k.presentation(), 0.5);
        check_pruning(nonzero_k.presentation(), 0.0);
        using Grid = gl::GridGradedSparseMatrix<double, 2, int>;
        check_pruning(Grid(nonzero_k.presentation()), 0.5);
        check_pruning(sd::convert_matrix_storage<gl::CSCStorage<int>>(nonzero_k.presentation()), 0.5);
        // A nonuniform physical grid exercises fractional shifts and reindexing.
        check_pruning(Grid(sd::Mat(0, 2, {}, {}, {{0.25, -2}, {1.25, -1}})), 0.5);
        // Every actual shifted map is zero, but the original lift quotient is
        // nonzero. Its convention still requires K = ker(eta), not K = 0.
        check_pruning(sd::Mat(2, 2, {{0}, {1}}, {{0.5, 0}, {1.5, 0}}, {{0, 0}, {1, 0}}), 1.0);
        check_pruning(sd::Mat(0, 0, {}, {}, {}), 0.5);
        if (argc == 3) check_pruning(sd::Module(argv[1]).presentation(), std::stod(argv[2]), false);
        std::cout << "PASS: all 1024 small-fixture switch combinations, pruning pairs, exact isomorphism, grids and CSC\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
