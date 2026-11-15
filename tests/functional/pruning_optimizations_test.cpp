// Independent switches preserve the pruning pair and output. The full-End
// control may visit extra passes when the quotient path returns immediately.
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
    sd::PruningOptions options;
    options.projected_preimages = mask & 1;
    options.amortize_reductions = mask & 2;
    options.full_end_space = mask & 4;
    options.composition_generators = (mask & 8) || (mask & 64);
    options.composition_budget.sparse_nilpotent = mask & 8;
    options.composition_budget.complement_orbit = mask & 64;
    options.adaptive_restriction = mask & 16;
    options.reuse_target_reductions = mask & 32;
    return options;
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
    const auto baseline_pair = sd::pruning_pair(*parent, epsilon, false, &baseline_profile, options_for(0));
    const auto baseline = sd::pruning(Module(matrix), epsilon, false, nullptr, options_for(0));
    // All three Hom algorithms must feed the same pruning pair, including
    // the composition-generator branch which previously forced exactness.
    for (auto algorithm : {gl::HomAlgorithm::optimised, gl::HomAlgorithm::full_restriction, gl::HomAlgorithm::hom_exactness}) {
        for (bool composition : {false, true}) {
            sd::PruningOptions options;
            options.hom_algorithm = algorithm;
            options.composition_generators = composition;
            const auto pair = sd::pruning_pair(*parent, epsilon, false, nullptr, options);
            require(pair.first.equals(baseline_pair.first) && pair.second.equals(baseline_pair.second),
                    "Hom algorithm changed the pruning pair");
        }
    }
    std::vector<unsigned> masks;
    if (exhaustive) {
        for (unsigned mask = 0; mask < 128; ++mask) masks.push_back(mask);
    } else {
        // Larger inputs cover full-End/composition/adaptive/reuse combinations with
        // two kernel/batching settings; small fixtures cover all 128 combinations.
        for (unsigned added = 0; added < 32; ++added)
            for (unsigned prior : {0u, 3u}) masks.push_back(prior | (added << 2));
    }
    for (unsigned mask : masks) {
        const auto options = options_for(mask);
        sd::PruningProfile profile;
        const auto pair = sd::pruning_pair(*parent, epsilon, false, &profile, options);
        require(pair.first.equals(baseline_pair.first), "I changed under an optimization");
        require(pair.second.equals(baseline_pair.second), "K changed under an optimization");
        if (!options.full_end_space)
            require(profile.i_iterations == baseline_profile.i_iterations &&
                    profile.k_iterations == baseline_profile.k_iterations, "Iteration counts changed");
        require(profile.i_restriction_kept <= profile.i_restriction_candidates &&
                profile.k_restriction_kept <= profile.k_restriction_candidates, "Adaptive restriction increased its family");
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
            auto reduced = sd::pruning_reduce_target(std::move(combined));
            require(reduced.row_degrees == map.row_degrees, "Target preparation changed the ambient basis");
            sd::Submodule reused(domain, sd::pruning_preimage_with_target(map, reduced, projected));
            require(reused.equals(expected), "Reused target changed the preimage");
        }
    }
}

void check_adaptive_restriction_grades() {
    auto parent = std::make_shared<const sd::Module>(sd::Mat(0, 1, {}, {}, {{0, 0}}));
    sd::Mat identity(1, 1, {{0}}, {{0, 0}}, {{0, 0}});
    std::vector<sd::Hom> maps{sd::Hom(parent, parent, identity), sd::Hom(parent, parent, identity)};
    sd::Mat target(1, 1, {{0}}, {{1, 1}}, {{0, 0}});
    require(sd::pruning_adaptive_restriction(maps, identity, target) == std::vector<std::size_t>{0},
            "A later target generator incorrectly killed an earlier map restriction");
    sd::Mat later(1, 1, {{0}}, {{2, 2}}, {{0, 0}});
    require(sd::pruning_adaptive_restriction(maps, later, target).empty(),
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
    for (unsigned added = 0; added < 32; ++added) for (unsigned prior : {0u, 3u}) {
        sd::PruningProfile profile;
        const auto options = options_for(prior | (added << 2));
        const auto result = sd::pruning_pair(*parent, 0.5, false, &profile, options);
        require(result.first.equals(expected.first) && result.second.equals(expected.second),
                "New switches changed pruning in another coordinate dimension");
        if (!options.full_end_space)
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
    auto free_parent = std::make_shared<const sd::Module>(free);
    for (bool full : {false, true}) for (bool composition : {false, true}) {
        sd::PruningOptions options;
        options.full_end_space = full;
        options.composition_generators = composition;
        sd::PruningProfile profile;
        auto [eta, i_maps, k_maps] = sd::prepare_pruning_pair(*free_parent, 0, &profile, options);
        require(profile.original_maps == (full ? 4 : 0),
                "Full-End control did not bypass the quotient");
        if (full && !composition)
            require(i_maps.size() == 4, "Full-End control did not retain the complete map basis");
    }
    // Exercise quotient coordinates beyond a machine word: 10 incomparable
    // free generators give dim Hom=100, dim End=10, hence dim Q=90.
    gl::vec<gl::r2degree> degrees;
    for (int j = 0; j < 10; ++j) degrees.push_back({double(j), double(9-j)});
    auto parent = std::make_shared<const sd::Module>(sd::Mat(0, 10, {}, {}, degrees));
    const auto eta = sd::Hom::canonical_shift(parent, {10, 10});
    sd::PruningProfile profile;
    const std::pair<std::vector<gl::Homomorphism<sd::Mat>>, std::vector<gl::Homomorphism<sd::Mat>>>
        families = sd::pruning_composition_generators(
            eta, gl::r2degree{10, 10}, gl::HomAlgorithm::hom_exactness, &profile);
    require(profile.original_maps == 90 && profile.quotient_dimension == 90,
            "Incorrect large Hom quotient dimension");
    require(profile.end_dimension == 10 && profile.i_maps == families.first.size() &&
            profile.k_maps == families.second.size(), "Composition statistics differ from map families");
    require(!families.first.empty() && families.first.size() < 90 &&
            !families.second.empty() && families.second.size() < 90,
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

void require_composition_span(const std::vector<sd::Hom>& basis, const std::vector<sd::Hom>& end,
                              const std::vector<sd::Hom>& maps, const gl::r2degree& delta,
                              const std::vector<sd::Hom>& denominator, bool left) {
    if (basis.empty()) { require(maps.empty(), "Empty quotient gained generators"); return; }
    const auto actions = gl::endomorphism_composition_action(basis, end, delta, left, denominator);
    const auto coordinates = gl::homomorphism_quotient_coordinates(basis, maps, denominator);
    gl::SparseMatrix<std::ptrdiff_t> generated(0, basis.size());
    for (const auto& action : actions) {
        const auto product = action * coordinates;
        for (std::ptrdiff_t c = 0; c < product.get_num_cols(); ++c)
            generated.append_col(product.get_col(c));
    }
    generated.compute_num_cols();
    generated.column_reduction_triangular(true);
    require(std::size_t(generated.get_num_cols()) == basis.size(),
            "Generators failed an independent complete-action spanning check");
}

// Independently build the complete action only on tiny fixtures, then check
// that every budgeted output generates ALL quotient coordinates on both sides.
void check_budgeted_composition_certificates() {
    const std::vector<sd::Mat> examples{
        sd::Mat(0, 3, {}, {}, {{0,0}, {1,1}, {2,2}}),
        sd::Mat(0, 4, {}, {}, {{0,1}, {0,1}, {1,0}, {1,0}}),
        sd::Mat(2, 3, {{0,1}, {1,2}}, {{2,2}, {3,3}}, {{0,0}, {1,1}, {2,2}}),
        // Nonminimal cover: a same-birth relation changes possible lifts.
        sd::Mat(1, 3, {{0,1}}, {{0,0}}, {{0,0}, {0,0}, {1,1}})
    };
    for (std::size_t example = 0; example < examples.size(); ++example) {
        auto parent = std::make_shared<const sd::Module>(examples[example]);
        const gl::r2degree delta{3,3};
        const auto eta = sd::Hom::canonical_shift(parent, delta);
        const auto end = gl::module_endomorphism_basis(parent, gl::HomAlgorithm::hom_exactness);
        const auto hom = gl::module_hom_space_basis(parent, eta.target(), gl::HomAlgorithm::hom_exactness);
        for (bool quotient : {false, true}) {
            std::vector<sd::Hom> denominator;
            if (quotient) for (const auto& e : end) denominator.push_back(e.compose(eta));
            const auto basis = gl::homomorphism_quotient_basis(hom, denominator);
            for (int setting = 0; setting < 5; ++setting) {
                gl::CompositionGeneratorOptions options;
                options.complement_orbit = false;
                if (setting == 1) options.max_products = 0;
                if (setting == 2) { options.max_products = 3; options.max_pair_visits = 5; }
                if (setting == 3) options.max_seeds = 0;
                if (setting == 4) options.max_product_terms = 0;
                std::array<gl::CompositionGeneratorStats, 2> stats;
                const auto result = gl::budgeted_composition_generators(basis, end, denominator, options, &stats);
                for (std::size_t side = 0; side < 2; ++side) {
                    const auto& maps = side ? result.second : result.first;
                    require(maps.size() <= basis.size(), "Budgeted search increased the family size");
                    require(stats[side].products <= options.max_products &&
                            stats[side].pair_visits <= options.max_pair_visits, "Composition budget exceeded");
                    if (setting == 1) require(maps.size() == basis.size(), "Zero budget lost basis vectors");
                    if (setting == 4) require(maps.size() == basis.size(), "Rejected products removed basis vectors");
                    for (std::size_t j = 0; j < stats[side].seeds; ++j) {
                        std::size_t entries = 0;
                        for (int c = 0; c < maps[j].generator_lift().get_num_cols(); ++c)
                            entries += maps[j].generator_lift().column(c).size();
                        require(entries <= options.max_seed_terms, "Retained candidate exceeded the density cap");
                    }
                    if (example == 1 && quotient && setting == 0)
                        require(maps.size() == 1, "Disjoint-support candidates missed the two identity cross blocks");
                    if (example == 1 && setting != 1)
                        require(stats[side].nilpotent_endomorphisms == 0,
                                "Equal-degree off-diagonal entries were mistaken for lowering maps");
                    if (example == 0 && setting == 3 && quotient)
                        require(stats[side].nilpotent_rank > 0 && maps.size() < basis.size(),
                                "Nilpotent actions failed to remove generators without retained seeds");
                    require_composition_span(basis, end, maps, delta, denominator, side == 0);
                }
            }
            // Independent switches, partial orbits, and chaining the original
            // coordinate certificate (the sparse output need not be a VS basis).
            const auto prepared = basis.empty() ? nullptr :
                std::make_unique<gl::PreparedCompositionCoordinates<sd::Mat>>(basis, denominator);
            for (unsigned method = 0; method < 4; ++method) for (int budget = 0; budget < 9; ++budget) {
                gl::CompositionGeneratorOptions options;
                options.sparse_nilpotent = method & 1;
                options.complement_orbit = method & 2;
                if (budget == 1) options.orbit_max_pair_visits = 1;
                if (budget == 2) options.orbit_max_seed_terms = 0;
                if (budget == 3) options.orbit_max_product_terms = 0;
                if (budget == 4) { options.orbit_rounds = 6; options.orbit_seed_mode = "sum"; }
                if (budget == 5) { options.orbit_rounds = 12; options.orbit_seed_mode = "random"; options.orbit_random_seed = 42; }
                if (budget == 6) options.orbit_rounds = 0;
                if (budget == 7) { options.orbit_min_gain = basis.size() + 1; options.orbit_patience = 2; }
                if (budget == 8) { options.orbit_max_product_terms = 0; options.orbit_patience = 0; options.orbit_rounds = 5; }
                std::array<gl::CompositionGeneratorStats, 2> stats;
                const auto result = gl::budgeted_composition_generators(
                    basis, end, denominator, options, &stats, prepared.get());
                if (budget == 5) {
                    const auto repeat = gl::budgeted_composition_generators(basis, end, denominator, options, nullptr, prepared.get());
                    for (std::size_t side = 0; side < 2; ++side) {
                        const auto& a = side ? result.second : result.first;
                        const auto& b = side ? repeat.second : repeat.first;
                        require(a.size() == b.size(), "Fixed random seed changed the generator count");
                        for (std::size_t j = 0; j < a.size(); ++j)
                            for (int c = 0; c < a[j].generator_lift().get_num_cols(); ++c)
                                require(a[j].generator_lift().get_col(c) == b[j].generator_lift().get_col(c),
                                        "Fixed random seed changed a generator");
                    }
                }
                for (std::size_t side = 0; side < 2; ++side) {
                    const auto& maps = side ? result.second : result.first;
                    require(maps.size() <= basis.size(), "Compression increased the family size");
                    const auto visit_limit = options.orbit_max_pair_visits ? options.orbit_max_pair_visits : options.orbit_rounds * end.size();
                    require(stats[side].orbit_pair_visits <= visit_limit, "Complement-orbit exceeded its visit budget");
                    require(stats[side].orbit_attempts <= options.orbit_rounds, "Orbit round limit exceeded");
                    if (budget == 7) require(stats[side].orbit_attempts <= 2, "Low-gain patience did not stop the search");
                    if (budget == 8 && method == 2 && !basis.empty() && !end.empty())
                        require(stats[side].orbit_attempts == 5,
                                "Disabled patience stopped rejected candidates early");
                    require(stats[side].largest_orbit_seed <= options.orbit_max_seed_terms,
                            "Complement-orbit exceeded its density budget");
                    if (method && !basis.empty()) {
                        require(maps.size() <= stats[side].sparse_generators,
                                "Adding the orbit stage increased the sparse output");
                        if (budget == 2 || budget == 3 || budget == 6 || budget == 8) require(maps.size() == stats[side].sparse_generators,
                                                "Skipped orbits changed the generating family");
                    }
                    if (method == 0) require(maps.size() == basis.size(), "Disabled methods lost basis maps");
                    require_composition_span(basis, end, maps, delta, denominator, side == 0);
                }
            }

        }
    }
}

void check_empty_hom_quotient() {
    // The first example had a nonempty lift quotient; the second did not.
    // Both have Q=0 and a zero shift map, hence I=M and K=M.
    for (const auto& matrix : {
            sd::Mat(2, 2, {{0}, {1}}, {{0.5,0}, {1.5,0}}, {{0,0}, {1,0}}),
            sd::Mat(1, 1, {{0}}, {{0.5,0}}, {{0,0}})}) {
        auto parent = std::make_shared<const sd::Module>(matrix);
        const auto eta = sd::Hom::canonical_shift(parent, {2,2});
        require(gl::End_d_0(eta).empty(), "Zero Hom quotient retained a lift");
        for (unsigned mask : {0u, 4u, 8u, 16u, 32u, 63u}) {
            const auto pair = sd::pruning_pair(*parent, 1.0, false, nullptr, options_for(mask));
            require(pair.first.equals(parent->whole_submodule()), "Empty quotient changed I");
            require(pair.second.equals(parent->whole_submodule()), "Empty quotient lost the shift kernel");
        }
        for (bool quick : {false, true}) {
            auto copy = matrix;
            const auto pair = sd::pruning_pair(copy, 1.0, quick);
            require(sd::Submodule(parent, pair.second).equals(parent->whole_submodule()),
                    "Matrix pruning lost the empty-quotient shift kernel");
            copy = matrix;
            sd::PruningProfile profile;
            const auto result = sd::pruning_profiled(copy, 1.0, quick, &profile);
            require(sd::Module(result).dimension_at({1,1}) == 0,
                    "Profiled pruning lost the empty-quotient shift kernel");
        }
    }
}

int main(int argc, char** argv) {
    try {
        check_projected_preimages();
        check_adaptive_restriction_grades();
        check_other_dimensions<1>();
        check_other_dimensions<3>();
        check_composition_generators();
        check_budgeted_composition_certificates();
        check_empty_hom_quotient();
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
        // nonzero. Empty Q must still give K = ker(eta), not K = 0.
        check_pruning(sd::Mat(2, 2, {{0}, {1}}, {{0.5, 0}, {1.5, 0}}, {{0, 0}, {1, 0}}), 1.0);
        check_pruning(sd::Mat(0, 0, {}, {}, {}), 0.5);
        if (argc == 3) check_pruning(sd::Module(argv[1]).presentation(), std::stod(argv[2]), false);
        std::cout << "PASS: all 128 small-fixture switch combinations, pruning pairs, exact isomorphism, grids and CSC\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
