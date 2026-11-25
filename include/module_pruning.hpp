/**
 * @file module_pruning.hpp
 * @brief Storage-independent pruning through the high-level Module API.
 *
 * Existing R2 entry points delegate here; coordinate and grid matrices use the same
 * iteration without accessing degree storage or matrix columns directly.
 */
#pragma once

#include "types.hpp"
#include "epsilon.hpp"
#include "progress.hpp"
#include "pruning_options.hpp"
#include "pruning_generators.hpp"
#include <grlina/hom_interface.hpp>
#include <grlina/storage/csc.hpp>
#include <grlina/matrices/grid.hpp>
#include <grlina/matrices/coord.hpp>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <grlina/matrices/matrix_algorithms.hpp>

namespace stable_decomposition {

// Grid matrices distinguish geometric shifts from stored index degrees.
template <typename Matrix>
auto diagonal_translation(const Matrix&, double amount) ->
    typename Matrix::degree_type {
    using Degree = typename Matrix::degree_type;
    using Scalar = typename Degree::value_type;
    if constexpr (std::is_integral_v<Scalar>) {
        if (std::trunc(amount) != amount ||
            static_cast<long double>(amount) < std::numeric_limits<Scalar>::lowest() ||
            static_cast<long double>(amount) > std::numeric_limits<Scalar>::max())
            throw std::invalid_argument(
                "Pruning shift cannot be represented by the grid scalar type");
    }
    Degree result{};
    result.fill(static_cast<Scalar>(amount));
    return result;
}

// Mutating pruning runs on vector columns; public CSC objects are packed only
// at the boundary, retaining their physical degrees and resolution maps.
template <typename Matrix>
inline constexpr bool editable_conversion_required =
    std::is_same_v<graded_linalg::matrix_storage_t<Matrix>,
                   graded_linalg::CSCStorage<graded_linalg::matrix_index_t<Matrix>>>;

template <class Storage, class Matrix>
auto convert_matrix_storage(const Matrix& matrix) {
    using Index = graded_linalg::matrix_index_t<Matrix>;
    using Result = graded_linalg::R2GradedSparseMatrix<Index, Storage>;
    graded_linalg::array<Index> columns;
    columns.reserve(matrix.get_num_cols());
    for (Index i = 0; i < matrix.get_num_cols(); ++i)
        columns.push_back(matrix.get_col(i));
    Result result(matrix.get_num_cols(), matrix.get_num_rows(), columns,
                  matrix.col_degrees, matrix.row_degrees);
    result.compatibly_sorted = matrix.compatibly_sorted;
    return result;
}

template <class Storage, class Matrix>
auto convert_module_storage(const graded_linalg::Module<Matrix>& module) {
    using ResultMatrix = decltype(convert_matrix_storage<Storage>(module.presentation()));
    auto convert_complex = [](const auto& complex) {
        graded_linalg::vec<ResultMatrix> maps;
        for (const auto& map : complex.differentials())
            maps.push_back(convert_matrix_storage<Storage>(map));
        return graded_linalg::ChainComplex<ResultMatrix>(std::move(maps));
    };
    auto result = graded_linalg::Module<ResultMatrix>::from_projective_resolution(
        convert_complex(module.proj_res_),
        module.res_complete_
            ? graded_linalg::ResolutionCompleteness::complete
            : graded_linalg::ResolutionCompleteness::truncated);
    result.inj_res_ = convert_complex(module.inj_res_);
    return result;
}

/** Return the canonical shift, the I map family, and the K map family.
 * Without composition selection, the I family is used for both loops and
 * the separate K family is left empty. Statistics are written only to profile.
 */
template <typename Matrix>
std::tuple<graded_linalg::Homomorphism<Matrix>,
           PruningMapFamily<Matrix>, PruningMapFamily<Matrix>>
prepare_pruning_pair(const graded_linalg::Module<Matrix>& X,
                     const typename Matrix::degree_type& epsilon,
                     PruningProfile* profile,
                     const PruningOptions& options) {
    const auto eps_vec = pruning_shift(epsilon);
    graded_linalg::Homomorphism<Matrix> eta_X =
        graded_linalg::Homomorphism<Matrix>::canonical_shift(X, eps_vec);

    std::pair<PruningCoefficients, PruningCoefficients> families;
    {
        PruningTimer timer(profile, &PruningProfile::basis);
        families = timed_with_progress(options.full_end_space ? "End_d" : "End_d/0", [&] {
            if (options.composition_enabled())
                return pruning_composition_coefficients(eta_X, eps_vec, options.hom_algorithm, profile,
                    options.full_end_space, options.composition_budget, options.single_solve_quotient);

            auto source = eta_X.domain()->presentation();
            source.compute_rows_forward();
            PruningCoefficients quotient_basis;
            if (options.full_end_space) {
                quotient_basis = graded_linalg::hom_space_coefficients(source, eta_X.target()->presentation(), options.hom_algorithm);
            } else if (options.single_solve_quotient) {
                quotient_basis = graded_linalg::End_d_0_coefficients(eta_X, options.hom_algorithm).quotient;
            } else {
                const auto end = graded_linalg::hom_space_coefficients(source, source, options.hom_algorithm);
                auto quotient = graded_linalg::End_d_0_coefficients(eta_X, end, options.hom_algorithm);
                quotient_basis = std::move(quotient.quotient);
            }
            if (profile) {
                profile->original_maps = quotient_basis.get_num_cols();
                profile->quotient_dimension = options.full_end_space ? 0 : quotient_basis.get_num_cols();
                profile->end_dimension = 0;
                profile->i_maps = profile->k_maps = quotient_basis.get_num_cols();
            }
            return std::make_pair(std::move(quotient_basis), PruningCoefficients(0, 0));
        });
    }

    auto maps = measure(profile, &PruningProfile::map_lifts, [&] { return std::make_pair(
        PruningMapFamily<Matrix>(std::move(families.first), eta_X.domain(), eta_X.target()),
        PruningMapFamily<Matrix>(std::move(families.second), eta_X.domain(), eta_X.target())); });
    return {std::move(eta_X), std::move(maps.first), std::move(maps.second)};
}

/** Shrink I by repeated intersections until a full pass leaves it unchanged. */
template <typename Matrix, class MapFamily>
std::pair<graded_linalg::Submodule<Matrix>, int>
compute_pruning_i(const graded_linalg::Homomorphism<Matrix>& eta_X,
                  const MapFamily& B_I,
                  PruningProfile* profile,
                  const PruningOptions& options, bool progress = true,
                  std::ostream* restriction_trace = nullptr, bool time_steps = false) {
    using Submodule = graded_linalg::Submodule<Matrix>;
    using Hom = graded_linalg::Homomorphism<Matrix>;
    if (options.threads < 1) throw std::invalid_argument("Thread count must be positive");
    if (options.parallel_i && options.update_i_target)
        throw std::invalid_argument("Parallel I requires a fixed target per pass");
    Submodule I = eta_X.domain()->whole_submodule();
    int iteration_I = 0;
    record_pruning_generators(profile, "I", "initial", 0, 0, I.num_of_gens());
    PruningTimer phase(profile, &PruningProfile::i_phase);
    auto shifted_image = [&] {
        auto image = measure(profile, &PruningProfile::i_shift, [&] {
            return eta_X.image(I, false);
        });
        if (options.minimize_i_target)
            measure(profile, &PruningProfile::i_target_reduce, [&] {
                image.minimize_generators(true, true);
            });
        return image;
    };
    for (;;) {
        std::size_t counter = 0;
        iteration_I++;
        if (profile)
            ++profile->i_iterations;
        const Submodule I_ = I;
        auto J = shifted_image();
        record_pruning_generators(profile, "J", options.minimize_i_target ? "reduced" : "raw",
                                  iteration_I, 0, J.num_of_gens());
        std::optional<Matrix> target_generators;
        auto prepare_target = [&] {
            target_generators = measure(profile, &PruningProfile::i_target_prepare, [&] {
                auto combined = eta_X.target()->presentation();
                combined.append_matrix(J.generator_map().generator_lift());
                if (options.reuse_target_reductions)
                    return pruning_reduce_target(std::move(combined));
                return combined;
            });
        };
        if ((options.adaptive_restriction || options.reuse_target_reductions) && !B_I.empty())
            prepare_target();
        std::vector<std::size_t> selected_I;
        PruningTiming restriction_timing;
        if (options.adaptive_restriction) {
            selected_I = measure(profile, &PruningProfile::i_adaptive_restriction, [&] {
                PruningTimer timer(time_steps ? &restriction_timing : nullptr);
                if (B_I.empty()) return std::vector<std::size_t>{};
                // Select using the pass-start I and J in both update modes.
                // Skipped constraints are reconsidered next pass; at convergence
                // I and J did not change, so this family certifies every constraint.
                return pruning_adaptive_restriction(B_I, I_.generator_map().generator_lift(), *target_generators, options.threads);
            });
            if (profile) {
                profile->i_restriction_candidates += B_I.size();
                profile->i_restriction_kept += selected_I.size();
                profile->i_restriction_max = std::max(profile->i_restriction_max, selected_I.size());
            }
        } else {
            selected_I.resize(B_I.size());
            std::iota(selected_I.begin(), selected_I.end(), std::size_t(0));
        }
        if (restriction_trace) {
            *restriction_trace << "restriction\tI\t" << iteration_I << '\t'
                               << B_I.size() << '\t' << selected_I.size();
            if (time_steps) *restriction_trace << '\t' << restriction_timing.seconds;
            *restriction_trace << '\n' << std::flush;
        }
        auto apply_constraint = [&](Submodule& current, std::size_t map_index, PruningProfile* timing) {
            const Hom& f = B_I[map_index];
            // If i generates I, I intersect f^-1(J) = i((f*i)^-1(J)).
            // Work on the free generator basis; no presentation of I is needed.
            auto restricted = measure(timing, &PruningProfile::i_image, [&] {
                return f.generator_lift() * current.generator_map().generator_lift();
            });
            auto inverse = measure(timing, &PruningProfile::i_preimage, [&] {
                if (options.reuse_target_reductions)
                    return pruning_preimage_with_target(restricted, *target_generators,
                        options.projected_preimages);
                return graded_linalg::inverse_image_generators(
                    restricted, f.target()->presentation(),
                    J.generator_map().generator_lift(), options.projected_preimages,
                    true);
            });
            auto next_I = measure(timing, &PruningProfile::i_image, [&] {
                return Submodule(current.parent(), current.generator_map().generator_lift() * inverse);
            });
            record_pruning_generators(timing, "I", "raw", iteration_I, map_index + 1, next_I.num_of_gens());
            measure(timing, &PruningProfile::i_reduce, [&] {
                next_I.minimize_generators(true, true);
            });
            const bool changed = options.update_i_target && !pruning_same_generators(
                current.generator_map().generator_lift(), next_I.generator_map().generator_lift());
            current = std::move(next_I);
            record_pruning_generators(timing, "I", "reduced", iteration_I, map_index + 1, current.num_of_gens());
            return changed;
        };
        if (options.parallel_i && options.threads > 1 && selected_I.size() > 1) {
            const auto workers = std::min(selected_I.size(), static_cast<std::size_t>(options.threads));
            std::vector<Submodule> partial(workers, I_);
            std::vector<PruningProfile> timings(workers);
            for (auto& timing : timings) if (profile) {
                timing.trace_step = profile->trace_step;
                timing.trace_generators = profile->trace_generators;
            }
            // Fixed J makes each group independent. Each worker shrinks its own
            // copy of pass-start I; intersection enforces every group's constraints.
            graded_linalg::parallel_for(workers, options.threads, [&](std::size_t worker) {
                const auto begin = selected_I.size() * worker / workers;
                const auto end = selected_I.size() * (worker + 1) / workers;
                for (auto j = begin; j < end; ++j)
                    apply_constraint(partial[worker], selected_I[j], profile ? &timings[worker] : nullptr);
            });
            // Merge statistics after joining: timer counters are never shared by workers.
            if (profile) for (const auto& timing : timings)
                for (auto field : {&PruningProfile::i_image, &PruningProfile::i_preimage, &PruningProfile::i_reduce}) {
                    (profile->*field).calls += (timing.*field).calls;
                    (profile->*field).seconds += (timing.*field).seconds;
                }
            I = std::move(partial.front());
            for (std::size_t worker = 1; worker < workers; ++worker) {
                if (pruning_same_generators(I.generator_map().generator_lift(),
                                           partial[worker].generator_map().generator_lift())) continue;
                I = measure(profile, &PruningProfile::i_intersection, [&] {
                    return I.intersection(partial[worker], true, options.projected_preimages, true);
                });
            }
        } else {
            bool target_changed = false;
            for (auto map_index : selected_I) {
                ++counter;
                if (progress && (counter == 1 ||
                    counter * 100 / selected_I.size() != (counter - 1) * 100 / selected_I.size()))
                    print_progress(iteration_I, counter, selected_I.size());
                if (options.update_i_target && target_changed) {
                    J = shifted_image();
                    record_pruning_generators(profile, "J", options.minimize_i_target ? "reduced" : "raw",
                                              iteration_I, map_index + 1, J.num_of_gens());
                    if (options.reuse_target_reductions) prepare_target();
                }
                target_changed = apply_constraint(I, map_index, profile);
            }
        }
        record_pruning_generators(profile, "I", "pass_end", iteration_I, 0, I.num_of_gens());
        if (progress) std::cout << std::endl;
        // I decreases, so only the reverse containment needs checking.
        if (measure(profile, &PruningProfile::i_equal, [&] {
                return I.contains(I_);
            }))
            break;
    }
    record_pruning_generators(profile, "I", "final", iteration_I, 0, I.num_of_gens());
    return {std::move(I), iteration_I};
}

/** Grow K by repeated sums and preimages until a full pass leaves it unchanged. */
template <typename Matrix, class MapFamily>
std::pair<graded_linalg::Submodule<Matrix>, int>
compute_pruning_k(const graded_linalg::Homomorphism<Matrix>& eta_X,
                  const MapFamily& B_K,
                  PruningProfile* profile,
                  const PruningOptions& options, bool progress = true,
                  std::ostream* restriction_trace = nullptr, bool time_steps = false) {
    using Submodule = graded_linalg::Submodule<Matrix>;
    using Hom = graded_linalg::Homomorphism<Matrix>;
    if (options.threads < 1) throw std::invalid_argument("Thread count must be positive");
    Submodule K = eta_X.domain()->zero_submodule();
    int iteration_K = 0;
    record_pruning_generators(profile, "K", "initial", 0, 0, K.num_of_gens());

    PruningTimer phase(profile, &PruningProfile::k_phase);
    for (;;) {
        std::size_t counter = 0;
        iteration_K++;
        if (profile)
            ++profile->k_iterations;
        const Submodule K_ = K;
        auto L = measure(profile, &PruningProfile::k_image, [&] {
            return eta_X.image(K, true, true);
        });
        record_pruning_generators(profile, "L", "reduced", iteration_K, 0, L.num_of_gens());
        std::vector<std::size_t> selected_K;
        PruningTiming restriction_timing;
        if (options.adaptive_restriction) {
            selected_K = measure(profile, &PruningProfile::k_adaptive_restriction, [&] {
                PruningTimer timer(time_steps ? &restriction_timing : nullptr);
                auto target = eta_X.target()->presentation();
                target.append_matrix(L.generator_map().generator_lift());
                return pruning_adaptive_restriction(B_K, K.generator_map().generator_lift(), target, options.threads);
            });
            if (profile) {
                profile->k_restriction_candidates += B_K.size();
                profile->k_restriction_kept += selected_K.size();
                profile->k_restriction_max = std::max(profile->k_restriction_max, selected_K.size());
            }
        } else {
            selected_K.resize(B_K.size());
            std::iota(selected_K.begin(), selected_K.end(), std::size_t(0));
        }
        if (restriction_trace) {
            *restriction_trace << "restriction\tK\t" << iteration_K << '\t'
                               << B_K.size() << '\t' << selected_K.size();
            if (time_steps) *restriction_trace << '\t' << restriction_timing.seconds;
            *restriction_trace << '\n' << std::flush;
        }
        // K and its initial image stay fixed. Each worker accumulates its own
        // image sum, keeping the existing bounded reductions within each group.
        auto accumulate_images = [&](Submodule& sum, std::size_t begin, std::size_t end,
                                     PruningProfile* timing, bool show_progress) {
            int pending = 0;
            auto accumulated = options.amortize_reductions
                ? sum.generator_map().generator_lift() : Matrix{};
            std::size_t current_map = 0;
            auto flush = [&] {
                sum = Submodule(sum.parent(), std::move(accumulated));
                measure(timing, &PruningProfile::k_reduce, [&] {
                    sum.minimize_generators(true, true);
                });
                record_pruning_generators(timing, "L", "reduced", iteration_K, current_map, sum.num_of_gens());
                accumulated = sum.generator_map().generator_lift();
                pending = 0;
            };
            for (auto j = begin; j < end; ++j) {
                const auto map_index = selected_K[j];
                current_map = map_index + 1;
                if (show_progress) {
                    ++counter;
                    if (counter == 1 || counter * 100 / selected_K.size() !=
                                        (counter - 1) * 100 / selected_K.size())
                        print_progress(iteration_K, counter, selected_K.size());
                }
                if (options.amortize_reductions && K.num_of_gens() == 0) continue;
                const Hom& f = B_K[map_index];
                auto image = measure(timing, &PruningProfile::k_image, [&] {
                    return f.image(K, false, true);
                });
                if (options.amortize_reductions) {
                    measure(timing, &PruningProfile::k_sum, [&] {
                        accumulated.append_matrix(image.generator_map().generator_lift());
                    });
                    record_pruning_generators(timing, "L", "raw", iteration_K, map_index + 1, accumulated.get_num_cols());
                    // ponytail: at most 16 images per worker between reductions;
                    // enlarge batches only if repeated reduction remains dominant.
                    if (++pending == 16) flush();
                } else {
                    sum = measure(timing, &PruningProfile::k_sum, [&] {
                        return sum.sum(image, true);
                    });
                    record_pruning_generators(timing, "L", "reduced", iteration_K, map_index + 1, sum.num_of_gens());
                }
            }
            if (options.amortize_reductions && pending) flush();
        };
        if (options.parallel_k && options.threads > 1 && selected_K.size() > 1 && K.num_of_gens() > 0) {
            const auto workers = std::min(selected_K.size(), static_cast<std::size_t>(options.threads));
            std::vector<Submodule> partial(workers, L);
            std::vector<PruningProfile> timings(workers);
            for (auto& timing : timings) if (profile) {
                timing.trace_step = profile->trace_step;
                timing.trace_generators = profile->trace_generators;
            }
            graded_linalg::parallel_for(workers, options.threads, [&](std::size_t worker) {
                accumulate_images(partial[worker], selected_K.size() * worker / workers,
                    selected_K.size() * (worker + 1) / workers,
                    profile ? &timings[worker] : nullptr, false);
            });
            if (profile) for (const auto& timing : timings)
                for (auto field : {&PruningProfile::k_image, &PruningProfile::k_sum, &PruningProfile::k_reduce}) {
                    (profile->*field).calls += (timing.*field).calls;
                    (profile->*field).seconds += (timing.*field).seconds;
                }
            L = std::move(partial.front());
            for (std::size_t worker = 1; worker < workers; ++worker) {
                if (pruning_same_generators(L.generator_map().generator_lift(),
                                           partial[worker].generator_map().generator_lift())) continue;
                // Fold reduced sums one at a time instead of storing one huge concatenation.
                L = measure(profile, &PruningProfile::k_merge, [&] {
                    return L.sum(partial[worker], true);
                });
                record_pruning_generators(profile, "L", "reduced", iteration_K, 0, L.num_of_gens());
            }
        } else {
            accumulate_images(L, 0, selected_K.size(), profile, progress);
        }
        // This implementation takes the preimage under eta_X rather than eta_I.
        // Consequently, K is represented as a submodule of X throughout the loop.
        K = measure(profile, &PruningProfile::k_preimage, [&] {
            auto result = eta_X.preimage(L, false, options.projected_preimages, true);
            record_pruning_generators(profile, "K", "raw", iteration_K, 0, result.num_of_gens());
            result.minimize_generators(true, true);
            record_pruning_generators(profile, "K", "reduced", iteration_K, 0, result.num_of_gens());
            return result;
        });
        record_pruning_generators(profile, "K", "pass_end", iteration_K, 0, K.num_of_gens());
        if (progress) std::cout << std::endl;
        // K increases: eta_X(K_) is included in L before taking its preimage.
        if (measure(profile, &PruningProfile::k_equal, [&] {
                return K_.contains(K);
            }))
            break;
    }
    record_pruning_generators(profile, "K", "final", iteration_K, 0, K.num_of_gens());
    return {std::move(K), iteration_K};
}

/** Coordinate and grid matrices share this high-level module algorithm.
 * profile optionally collects timing and operation counts; nullptr disables it.
 */
template <typename Matrix>
std::pair<graded_linalg::Submodule<Matrix>, graded_linalg::Submodule<Matrix>>
pruning_pair(const graded_linalg::Module<Matrix>& X,
             const typename Matrix::degree_type& epsilon,
             bool quick = false,
             PruningProfile* profile = nullptr,
             const PruningOptions& options = {}) {

    if (options.threads < 1) throw std::invalid_argument("Thread count must be positive");

    if constexpr (editable_conversion_required<Matrix>) {
        auto original = X.shared_from_this();
        using EditableStorage = graded_linalg::vec<graded_linalg::vec<graded_linalg::matrix_index_t<Matrix>>>;
        auto editable = convert_module_storage<EditableStorage>(X);
        auto parent = std::make_shared<const decltype(editable)>(std::move(editable));
        auto [I, K] =
            stable_decomposition::pruning_pair(*parent, epsilon, quick, profile, options);
        // Preserve the caller's exact parent identity. Pair results represent
        // embeddings; their optional presentation cache is computed on demand.
        return {graded_linalg::Submodule<Matrix>(
                    original,
                    convert_matrix_storage<graded_linalg::matrix_storage_t<Matrix>>(
                        I.generator_map().generator_lift())),
                graded_linalg::Submodule<Matrix>(
                    original,
                    convert_matrix_storage<graded_linalg::matrix_storage_t<Matrix>>(
                        K.generator_map().generator_lift()))};
    } else {

        (void)quick; // the high-level algorithm retains the existing ordinary iteration
        auto preparation = prepare_pruning_pair(X, epsilon, profile, options);
        const auto& eta_X = std::get<0>(preparation);
        const auto& B_I = std::get<1>(preparation);
        const auto& composition_k_maps = std::get<2>(preparation);
        const auto& B_K =
            options.composition_enabled() ? composition_k_maps : B_I;

        // Q = 0 imposes no I constraints, but K must still contain ker(eta).
        if (B_I.empty() && B_K.empty()) {
            auto I = X.whole_submodule();
            record_pruning_generators(profile, "I", "initial", 0, 0, I.num_of_gens());
            record_pruning_generators(profile, "K", "initial", 0, 0, 0);
            PruningTimer timer(profile, &PruningProfile::k_preimage);
            auto K = eta_X.kernel(false);
            record_pruning_generators(profile, "K", "raw", 0, 0, K.num_of_gens());
            K.minimize_generators();
            record_pruning_generators(profile, "K", "reduced", 0, 0, K.num_of_gens());
            record_pruning_generators(profile, "I", "final", 0, 0, I.num_of_gens());
            record_pruning_generators(profile, "K", "final", 0, 0, K.num_of_gens());
            return {std::move(I), std::move(K)};
        }

        auto phases = [&] {
            if (options.threads > 1) {
                auto i_options = options, k_options = options;
                // ponytail: fixed split; rebalance only if phase imbalance dominates.
                k_options.threads = options.threads / 2;
                i_options.threads -= k_options.threads;
                // Read-only parents/maps; I and K write disjoint profile fields.
                auto i = std::async(std::launch::async, [&] {
                    return compute_pruning_i(eta_X, B_I, profile, i_options, false);
                });
                auto k = compute_pruning_k(eta_X, B_K, profile, k_options, false);
                return std::make_pair(i.get(), std::move(k));
            }
            auto i = compute_pruning_i(eta_X, B_I, profile, options);
            auto k = compute_pruning_k(eta_X, B_K, profile, options);
            return std::make_pair(std::move(i), std::move(k));
        }();
        auto [I, iteration_I] = std::move(phases.first);
        auto [K, iteration_K] = std::move(phases.second);

        std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K"
                  << std::endl;
        // The standard estimate counts fixed-target passes; a fast pass can
        // propagate several constraints and does not give the same bound.
        if (!options.update_i_target) {
            auto bound = epsilon;
            for (auto& component : bound) component *= std::max(iteration_I, iteration_K);
            std::cout << "Output is " << format_epsilon(bound)
                      << "-interleaved with the input" << std::endl;
        }
        return {std::move(I), std::move(K)};
    }
}

template <typename Matrix>
graded_linalg::Module<Matrix> pruning(graded_linalg::Module<Matrix> X,
                                      const typename Matrix::degree_type& epsilon,
                                      bool quick = false,
                                      PruningProfile* profile = nullptr,
                                      const PruningOptions& options = {}) {
    if constexpr (editable_conversion_required<Matrix>) {
        using EditableStorage = graded_linalg::vec<graded_linalg::vec<graded_linalg::matrix_index_t<Matrix>>>;
        auto editable = convert_module_storage<EditableStorage>(X);
        auto result = stable_decomposition::pruning(
            std::move(editable), epsilon, quick, profile, options);
        return convert_module_storage<graded_linalg::matrix_storage_t<Matrix>>(result);
    } else {
        using Module = graded_linalg::Module<Matrix>;
        if (profile) {
            auto trace = std::move(profile->trace_step);
            auto generators = std::move(profile->trace_generators);
            *profile = {};
            profile->trace_step = std::move(trace);
            profile->trace_generators = std::move(generators);
        }
        PruningTimer total(profile, &PruningProfile::total);
        auto parent = std::make_shared<Module>(std::move(X));
        auto [I, K] = stable_decomposition::pruning_pair<Matrix>(*parent, epsilon, quick, profile, options);
        auto I_qK = measure(profile, &PruningProfile::quotient, [&I = I, &K = K] {
            return I.submodule_quotient(K);
        });
        measure(profile, &PruningProfile::presentation, [&] {
            I_qK.compute_presentation();
        });
        Module result = std::move(I_qK);
        measure(profile, &PruningProfile::minimize, [&] {
            result.minimize();
        });
        // Translate physical coordinates; grid indices alone do not encode epsilon.
        auto translation = epsilon;
        for (auto& component : translation) component = -component;
        result.shift(translation);
        return result;
    }
}

// Scalar overloads retain the diagonal API in any number of parameters.
template <typename Matrix>
auto prepare_pruning_pair(const graded_linalg::Module<Matrix>& X, double epsilon,
                          PruningProfile* profile, const PruningOptions& options) {
    (void)pruning_shift(epsilon);
    return prepare_pruning_pair(X, diagonal_translation(X.presentation(), epsilon), profile, options);
}
template <typename Matrix>
auto pruning_pair(const graded_linalg::Module<Matrix>& X, double epsilon, bool quick = false,
                  PruningProfile* profile = nullptr, const PruningOptions& options = {}) {
    (void)pruning_shift(epsilon);
    return stable_decomposition::pruning_pair(X, diagonal_translation(X.presentation(), epsilon), quick, profile, options);
}
template <typename Matrix>
auto pruning(graded_linalg::Module<Matrix> X, double epsilon, bool quick = false,
             PruningProfile* profile = nullptr, const PruningOptions& options = {}) {
    (void)pruning_shift(epsilon);
    const auto amount = diagonal_translation(X.presentation(), epsilon);
    return stable_decomposition::pruning(std::move(X), amount, quick, profile, options);
}
} // namespace stable_decomposition
