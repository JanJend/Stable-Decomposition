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
    options.threads = 1;
    options.projected_preimages = mask & 1;
    options.amortize_reductions = mask & 2;
    options.full_end_space = mask & 4;
    options.composition_generators = (mask & 8) || (mask & 64);
    options.composition_budget.nilpotent_reduction = mask & 8;
    options.composition_budget.complement_orbit = mask & 64;
    options.adaptive_restriction = mask & 16;
    options.reuse_target_reductions = mask & 32;
    options.single_solve_quotient = mask & 128;
    options.update_i_target = mask & 256;
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
    sd::PruningProfile serial_profile;
    sd::PruningOptions serial_options;
    serial_options.threads = 1;
    const auto serial_pair = sd::pruning_pair(*parent, epsilon, false, &serial_profile, serial_options);
    for (int threads : {sd::PruningOptions{}.threads, 2, 3, 4}) {
        sd::PruningOptions options;
        options.threads = threads;
        sd::PruningProfile profile;
        const auto pair = sd::pruning_pair(*parent, epsilon, false, &profile, options);
        require(pair.first.equals(serial_pair.first) && pair.second.equals(serial_pair.second),
                "Parallel pruning changed I or K");
        require(pair.first.parent() == parent && pair.second.parent() == parent,
                "Parallel pruning changed parent identity");
        require(profile.i_iterations == serial_profile.i_iterations &&
                profile.k_iterations == serial_profile.k_iterations &&
                profile.i_restriction_candidates == serial_profile.i_restriction_candidates &&
                profile.k_restriction_candidates == serial_profile.k_restriction_candidates &&
                profile.i_restriction_kept == serial_profile.i_restriction_kept &&
                profile.k_restriction_kept == serial_profile.k_restriction_kept &&
                profile.i_preimage.calls == serial_profile.i_preimage.calls &&
                profile.k_preimage.calls == serial_profile.k_preimage.calls,
                "Parallel pruning changed iteration counts, restrictions or profiling counts");
        const auto result = sd::pruning(Module(matrix), epsilon, false, nullptr, options);
        require(gl::is_isomorphic(coordinate_matrix(baseline.presentation()),
                                 coordinate_matrix(result.presentation())),
                "Parallel pruning changed the output module");
        options.parallel_i = true;
        for (bool minimize_target : {false, true}) {
            options.minimize_i_target = minimize_target;
            const auto parallel_pair = sd::pruning_pair(*parent, epsilon, false, nullptr, options);
            require(parallel_pair.first.equals(serial_pair.first) && parallel_pair.second.equals(serial_pair.second),
                    "Parallel fixed-target I or J minimization changed the pruning pair");
        }
    }
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
        for (unsigned mask = 0; mask < 512; ++mask) masks.push_back(mask);
    } else {
        // Larger inputs cover full-End/composition/adaptive/reuse combinations with
        // two kernel/batching settings; small fixtures cover all 512 combinations.
        for (unsigned added = 0; added < 128; ++added)
            for (unsigned prior : {0u, 3u}) masks.push_back(prior | (added << 2));
    }
    for (unsigned mask : masks) {
        auto options = options_for(mask);
        options.threads = 1 + mask % 4;
        sd::PruningProfile profile;
        const auto pair = sd::pruning_pair(*parent, epsilon, false, &profile, options);
        require(pair.first.equals(baseline_pair.first), "I changed under an optimization");
        require(pair.second.equals(baseline_pair.second), "K changed under an optimization");
        if (!options.full_end_space) {
            require(profile.k_iterations == baseline_profile.k_iterations, "K iteration counts changed");
            if (!options.update_i_target)
                require(profile.i_iterations == baseline_profile.i_iterations, "Standard I iteration counts changed");
        }
        require(profile.i_restriction_kept <= profile.i_restriction_candidates &&
                profile.k_restriction_kept <= profile.k_restriction_candidates, "Adaptive restriction increased its family");
        const auto result = sd::pruning(Module(matrix), epsilon, false, nullptr, options);
        require(result.presentation().is_graded_matrix(), "Ungraded optimized result");
        require(gl::is_isomorphic(coordinate_matrix(baseline.presentation()),
                                   coordinate_matrix(result.presentation())), "Pruning output changed");
    }
}

void check_i_target_updates() {
    // f sends e1 to e2, then g sends e0 to e1. With shift (2,2), the
    // updated target propagates both constraints in one pass instead of two.
    const gl::vec<gl::r2degree> births{{0,0}, {1,1}, {2,2}};
    auto parent = std::make_shared<const sd::Module>(sd::Mat(0, 3, {}, {}, births));
    const auto eta = sd::Hom::canonical_shift(parent, {2,2});
    const auto& shifted = eta.target()->presentation().row_degrees;
    const std::vector<sd::Hom> maps{
        sd::Hom(parent, eta.target(), sd::Mat(3, 3, {{}, {2}, {}}, births, shifted)),
        sd::Hom(parent, eta.target(), sd::Mat(3, 3, {{1}, {}, {}}, births, shifted))};
    const sd::Submodule expected(parent,
        sd::Mat(3, 3, {{0}, {1}, {2}}, {{2,2}, {2,2}, {2,2}}, births));
    for (bool adaptive : {false, true}) for (bool reuse : {false, true})
        for (bool projected : {false, true}) for (bool updated : {false, true})
        for (bool minimize_target : {false, true}) {
            auto options = options_for(0);
            options.adaptive_restriction = adaptive;
            options.reuse_target_reductions = reuse;
            options.projected_preimages = projected;
            options.update_i_target = updated;
            options.minimize_i_target = minimize_target;
            const auto result = sd::compute_pruning_i(eta, maps, nullptr, options, false);
            require(result.first.equals(expected), "I target mode changed the stable submodule");
            require(result.second == (updated ? 2 : 3), "I target mode did not change pass propagation");
            if (!updated) for (int threads : {1, 2, 4}) {
                options.parallel_i = true;
                options.threads = threads;
                sd::PruningProfile timing;
                const auto parallel = sd::compute_pruning_i(eta, maps, &timing, options, false);
                require(parallel.first.equals(expected) && parallel.second == result.second,
                        "Independent constraint groups changed I or its fixed-target passes");
                require(timing.i_preimage.calls == timing.i_reduce.calls && timing.i_preimage.calls > 0,
                        "Parallel I lost worker timings");
                require((timing.i_target_reduce.calls > 0) == minimize_target,
                        "J was minimized despite the requested setting");
                if (threads > 1 && !adaptive)
                    require(timing.i_intersection.calls > 0, "Parallel I did not merge different worker submodules");
            }
        }
    // Identity constraints leave I fixed. Fast I must reuse its image/target
    // while the changing-constraint case above still refreshes when necessary.
    auto free = std::make_shared<const sd::Module>(sd::Mat(0, 2, {}, {}, {{0,0}, {1,1}}));
    const auto shift = sd::Hom::canonical_shift(free, {2,2});
    std::vector<sd::Hom> identities(70, shift);
    auto options = options_for(0);
    options.update_i_target = options.reuse_target_reductions = true;
    sd::PruningProfile profile;
    const auto fixed = sd::compute_pruning_i(shift, identities, &profile, options, false);
    require(fixed.second == 1 && fixed.first.equals(free->whole_submodule()), "Identity constraints changed I");
    require(profile.i_image.calls == 2 * identities.size() && profile.i_shift.calls == 1 &&
            profile.i_target_prepare.calls == 1 && profile.i_target_reduce.calls == 0,
            "Fast I refreshed an unchanged image or target");
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
    sd::Mat later(1, 1, {{0}}, {{2, 2}}, {{0, 0}});
    for (int threads : {1, 2, 3, 4}) {
        require(sd::pruning_adaptive_restriction(maps, identity, target, threads) == std::vector<std::size_t>{0},
                "A later target generator incorrectly killed an earlier map restriction");
        require(sd::pruning_adaptive_restriction(maps, later, target, threads).empty(),
                "A target generator did not kill a restriction at a later grade");
    }
    for (std::size_t failing_job : {0u, 1u}) {
        bool caught = false;
        try {
            gl::parallel_for(4, 2, [&](std::size_t i) {
                if (i == failing_job) throw std::runtime_error("worker failure");
            });
        } catch (const std::runtime_error& error) {
            caught = std::string(error.what()) == "worker failure";
        }
        require(caught, "Parallel worker exception was lost");
    }
}

void check_coefficient_families_and_batches() {
    // Duplicates and F2 cancellations straddle several 32-map batch boundaries.
    auto parent = std::make_shared<const sd::Module>(sd::Mat(0, 3, {}, {}, {{0,0}, {1,0}, {0,1}}));
    const auto eta = sd::Hom::canonical_shift(parent, {3,3});
    sd::PruningCoefficients coefficients(97, 9);
    coefficients.resize_data(97);
    std::mt19937 random(572934);
    for (int c = 0; c < 97; ++c) {
        gl::vec<std::ptrdiff_t> entries;
        if (c % 3) for (int r = 0; r < 9; ++r) if (random() % 2) entries.push_back(r);
        coefficients.set_col(c, entries);
    }
    const auto maps = gl::homomorphisms_from_coefficients(coefficients, eta.domain(), eta.target());
    const sd::PruningMapFamily<sd::Mat> family(coefficients, eta.domain(), eta.target());
    const auto source = parent->whole_submodule().generator_map().generator_lift();
    sd::Mat target(2, 3, {{0,1}, {1,2}}, {{1,0}, {0,1}}, eta.target()->presentation().row_degrees);
    std::vector<sd::Mat> restrictions;
    for (std::size_t i = 0; i < maps.size(); ++i) {
        const auto lift = family[i].generator_lift();
        require(sd::pruning_same_generators(lift, maps[i].generator_lift()), "Coefficient family changed a lift");
        require(lift.row_degrees == maps[i].generator_lift().row_degrees, "Coefficient family changed target grades");
        restrictions.push_back(maps[i].generator_lift() * source);
    }
    const auto expected = gl::independent_matrix_family_indices_modulo_target_relations(restrictions, target);
    for (int threads : {1, 2, 4}) {
        require(sd::pruning_adaptive_restriction(family, source, target, threads) == expected,
                "Batched coefficient restrictions changed the independent indices");
        require(sd::pruning_adaptive_restriction(maps, source, target, threads) == expected,
                "Batched materialized restrictions changed the independent indices");
    }
    using Grid = gl::GridGradedSparseMatrix<double, 2, int>;
    auto grid_parent = std::make_shared<const gl::Module<Grid>>(Grid(parent->presentation()));
    const auto grid_eta = gl::Homomorphism<Grid>::canonical_shift(grid_parent, {3,3});
    const sd::PruningMapFamily<Grid> grid_family(coefficients, grid_eta.domain(), grid_eta.target());
    const auto grid_source = grid_parent->whole_submodule().generator_map().generator_lift();
    require(sd::pruning_adaptive_restriction(grid_family, grid_source, Grid(target), 2) == expected,
            "Batched restrictions changed physical grid degrees");
    auto later = source;
    later.col_degrees[0] = {0.5,0.5};
    require(!sd::pruning_same_generators(source, later), "Cache equality ignored generator degrees");
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
    for (unsigned added = 0; added < 64; ++added) for (unsigned prior : {0u, 3u}) {
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
        sd::Mat(1, 3, {{0,1}}, {{0,0}}, {{0,0}, {0,0}, {1,1}}),
        // No generator ever vanishes and every generator occurs in a relation:
        // the inexpensive ungraded check must reject all elementary factors.
        sd::Mat(2, 3, {{0,1}, {1,2}}, {{3,3}, {4,4}}, {{0,0}, {1,1}, {2,2}})
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
            // Check the production quotient representatives against the
            // independently computed full Hom space and canonical-shift image.
            for (auto algorithm : {gl::HomAlgorithm::optimised, gl::HomAlgorithm::full_restriction,
                                   gl::HomAlgorithm::hom_exactness}) for (bool single : {false, true}) {
                const auto families = sd::pruning_composition_generators(
                    eta, delta, algorithm, nullptr, !quotient, {}, single);
                require_composition_span(basis, end, families.first, delta, denominator, true);
                require_composition_span(basis, end, families.second, delta, denominator, false);
            }
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
                    require(stats[side].products + stats[side].elementary_certificates <= options.max_products &&
                            stats[side].certificate_visits <= options.max_certificate_visits,
                            "Elementary certificate budget exceeded");
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
                    if (example == 4)
                        require(stats[side].certificate_visits == 0 && stats[side].elementary_certificates == 0,
                                "Impossible elementary factors were searched");
                    require_composition_span(basis, end, maps, delta, denominator, side == 0);
                }
            }
            // Independent switches, partial orbits, and chaining the original
            // coordinate certificate (the nilpotent-reduction output need not be a VS basis).
            const auto prepared = basis.empty() ? nullptr :
                std::make_unique<gl::PreparedCompositionCoordinates<sd::Mat>>(basis, denominator);
            for (unsigned method = 0; method < 4; ++method) for (int budget = 0; budget < 9; ++budget) {
                gl::CompositionGeneratorOptions options;
                options.nilpotent_reduction = method & 1;
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
                        require(maps.size() <= stats[side].nilpotent_generators,
                                "Adding the orbit stage increased the nilpotent-reduction output");
                        if (budget == 2 || budget == 3 || budget == 6 || budget == 8) require(maps.size() == stats[side].nilpotent_generators,
                                                "Skipped orbits changed the generating family");
                    }
                    if (method == 0) require(maps.size() == basis.size(), "Disabled methods lost basis maps");
                    require_composition_span(basis, end, maps, delta, denominator, side == 0);
                }
            }

        }
    }
}

void check_elementary_composition_certificates() {
    const std::vector<gl::r2degree> births{{0,2}, {1,3}, {3,0}};
    const gl::r2degree delta{3,3};
    const std::vector<sd::Mat> examples{
        sd::Mat(0, 3, {}, {}, births),
        // E_01 violates the relation g_0 + g_1, although the row match exists.
        sd::Mat(1, 3, {{0,1}}, {{4,4}}, births),
        // E_01 is valid, but moving E_02 to row 1 violates the death of g_2.
        sd::Mat(3, 3, {{2}, {0}, {1}}, {{4,0}, {7,3}, {8,4}}, births)
    };
    for (std::size_t example = 0; example < examples.size(); ++example) {
        auto parent = std::make_shared<const sd::Module>(examples[example]);
        const auto eta = sd::Hom::canonical_shift(parent, delta);
        auto end = gl::module_endomorphism_basis(parent, gl::HomAlgorithm::hom_exactness);
        const auto hom = gl::module_hom_space_basis(parent, eta.target(), gl::HomAlgorithm::hom_exactness);
        std::vector<sd::Hom> denominator;
        for (const auto& e : end) denominator.push_back(e.compose(eta));
        const auto basis = gl::homomorphism_quotient_basis(hom, denominator);
        const auto target_births = gl::geometric_row_degrees(eta.target()->presentation());
        const sd::Mat t(3, 3, {{}, {0}, {}}, births, births);
        const sd::Mat h(3, 3, {{}, {}, {1}}, births, target_births);
        const sd::Mat f(3, 3, {{}, {}, {0}}, births, target_births);
        require(gl::is_homomorphism(examples[example], eta.target()->presentation(), f),
                "Elementary certificate fixture has an invalid F");
        require(gl::is_homomorphism(examples[example], examples[example], t) == (example != 1),
                "Elementary certificate fixture has the wrong T compatibility");
        require(gl::is_homomorphism(examples[example], eta.target()->presentation(), h) == (example != 2),
                "Elementary certificate fixture has the wrong H compatibility");
        if (example == 0) {
            // These one-dimensional Hom subspaces are invariant, but their
            // row/column-moved factors lie outside the supplied subspace.
            gl::CompositionGeneratorOptions isolated;
            isolated.complement_orbit = false;
            isolated.nilpotent_rounds = isolated.max_seeds = 0;
            const sd::Mat right_f(3, 3, {{}, {2}, {}}, births, target_births);
            for (const auto& lift : {f, right_f}) {
                const std::vector<sd::Hom> subspace{sd::Hom(parent, eta.target(), lift)};
                std::array<gl::CompositionGeneratorStats, 2> stats;
                const auto result = gl::budgeted_composition_generators(subspace, end, {}, isolated, &stats);
                require(result.first.size() == 1 && result.second.size() == 1 &&
                        stats[0].elementary_certificates == 0 && stats[1].elementary_certificates == 0,
                        "A factor outside the supplied Hom subspace removed a generator");
            }
            // Hide E_01 in a basis combination. The new pass must validate this
            // factor on demand instead of relying on a lowering End basis vector.
            for (auto& e : end) {
                auto& lift = e.generator_lift();
                if (lift.column(0).empty() && lift.get_col(1) == gl::vec<int>{0} && lift.column(2).empty())
                    lift.set_col(0, {0});
            }
        }
        for (std::size_t visits : {std::size_t(0), std::size_t(1), std::size_t(8192)}) {
            gl::CompositionGeneratorOptions options;
            options.complement_orbit = false;
            options.nilpotent_rounds = options.max_seeds = 0;
            options.max_certificate_visits = visits;
            std::array<gl::CompositionGeneratorStats, 2> stats;
            const auto result = gl::budgeted_composition_generators(basis, end, denominator, options, &stats);
            for (std::size_t side = 0; side < 2; ++side) {
                const auto& maps = side ? result.second : result.first;
                require(stats[side].products == 0 && stats[side].certificate_visits <= visits,
                        "Direct certificate performed a product or exceeded its visit cap");
                if (!visits) require(maps.size() == basis.size(), "Disabled certificates removed maps");
                if (example == 0 && visits == 8192)
                    require(stats[side].nilpotent_endomorphisms == 0 &&
                            stats[side].elementary_certificates > 0 && maps.size() < basis.size(),
                            "Direct certificate missed a factor hidden in End basis combinations");
                if (example == 1 || (example == 2 && side == 0))
                    require(stats[side].elementary_certificates == 0,
                            "Invalid elementary factor passed the relation checks");
                require_composition_span(basis, end, maps, delta, denominator, side == 0);
            }
        }
        if (example == 0) {
            gl::CompositionGeneratorOptions options;
            options.complement_orbit = false;
            options.max_products = 2; // One direct certificate fills the nilpotent half.
            std::array<gl::CompositionGeneratorStats, 2> stats;
            const auto result = gl::budgeted_composition_generators(basis, end, denominator, options, &stats);
            for (std::size_t side = 0; side < 2; ++side) {
                require(stats[side].elementary_certificates == 1 && stats[side].products == 0 &&
                        stats[side].coordinate_preparations == 0 && stats[side].seeds == 0,
                        "Filled direct-certificate budget still prepared coordinates or ordinary seeds");
                require_composition_span(basis, end, side ? result.second : result.first,
                                         delta, denominator, side == 0);
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
        check_i_target_updates();
        check_adaptive_restriction_grades();
        check_coefficient_families_and_batches();
        check_other_dimensions<1>();
        check_other_dimensions<3>();
        check_composition_generators();
        check_budgeted_composition_certificates();
        check_elementary_composition_certificates();
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
        std::cout << "PASS: all 512 small-fixture switch combinations, both I target modes, 1-4 threads, pruning pairs, exact isomorphism, grids and CSC\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
