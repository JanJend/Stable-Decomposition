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
           std::vector<graded_linalg::Homomorphism<Matrix>>,
           std::vector<graded_linalg::Homomorphism<Matrix>>>
prepare_pruning_pair(const graded_linalg::Module<Matrix>& X,
                     double epsilon,
                     PruningProfile* profile,
                     const PruningOptions& options) {
    (void)pruning_shift(epsilon); // shared finite/nonnegative validation
    const typename Matrix::degree_type eps_vec =
        diagonal_translation(X.presentation(), 2 * epsilon);
    graded_linalg::Homomorphism<Matrix> eta_X =
        graded_linalg::Homomorphism<Matrix>::canonical_shift(X, eps_vec);

    std::pair<std::vector<graded_linalg::Homomorphism<Matrix>>,
              std::vector<graded_linalg::Homomorphism<Matrix>>> families;
    {
        PruningTimer timer(profile ? &profile->basis : nullptr);
        families = timed_with_progress(options.full_end_space ? "End_d" : "End_d/0", [&] {
            if (options.composition_enabled())
                return pruning_composition_generators(eta_X, eps_vec, options.hom_algorithm, profile, options.full_end_space, options.composition_budget);

            std::vector<graded_linalg::Homomorphism<Matrix>> quotient_basis = options.full_end_space
                ? graded_linalg::module_hom_space_basis(eta_X.domain(), eta_X.target(), options.hom_algorithm)
                : graded_linalg::End_d_0(eta_X, true, options.hom_algorithm);
            if (profile) {
                profile->original_maps = quotient_basis.size();
                profile->quotient_dimension = options.full_end_space ? 0 : quotient_basis.size();
                profile->end_dimension = 0;
                profile->i_maps = profile->k_maps = quotient_basis.size();
            }
            return std::make_pair(std::move(quotient_basis),
                std::vector<graded_linalg::Homomorphism<Matrix>>{});
        });
    }

    return {std::move(eta_X), std::move(families.first), std::move(families.second)};
}

/** Shrink I by repeated intersections until a full pass leaves it unchanged. */
template <typename Matrix>
std::pair<graded_linalg::Submodule<Matrix>, int>
compute_pruning_i(const graded_linalg::Homomorphism<Matrix>& eta_X,
                  const std::vector<graded_linalg::Homomorphism<Matrix>>& B_I,
                  PruningProfile* profile,
                  const PruningOptions& options) {
    using Submodule = graded_linalg::Submodule<Matrix>;
    using Hom = graded_linalg::Homomorphism<Matrix>;
    Submodule I = eta_X.domain()->whole_submodule();
    int iteration_I = 0;
    PruningTimer phase(profile ? &profile->i_phase : nullptr);
    for (;;) {
        std::size_t counter = 0;
        iteration_I++;
        if (profile)
            ++profile->i_iterations;
        const Submodule I_ = I;
        auto J = measure(profile, &PruningProfile::i_image, [&] {
            return eta_X.image(I, true, true);
        });
        std::optional<Matrix> target_generators;
        if ((options.adaptive_restriction || options.reuse_target_reductions) && !B_I.empty()) {
            target_generators = measure(profile, &PruningProfile::i_target_prepare, [&] {
                auto combined = eta_X.target()->presentation();
                combined.append_matrix(J.generator_map().generator_lift());
                if (options.reuse_target_reductions)
                    return pruning_reduce_target(std::move(combined));
                return combined;
            });
        }
        std::vector<std::size_t> selected_I;
        if (options.adaptive_restriction) {
            selected_I = measure(profile, &PruningProfile::i_adaptive_restriction, [&] {
                if (B_I.empty()) return std::vector<std::size_t>{};
                // Use the pass-start I; J remains fixed through the pass.
                return pruning_adaptive_restriction(B_I, I_.generator_map().generator_lift(), *target_generators);
            });
            if (profile) {
                profile->i_restriction_candidates += B_I.size();
                profile->i_restriction_kept += selected_I.size();
            }
        } else {
            selected_I.resize(B_I.size());
            std::iota(selected_I.begin(), selected_I.end(), std::size_t(0));
        }
        for (auto map_index : selected_I) {
            const Hom& f = B_I[map_index];
            ++counter;
            // At most one update per percentage, with no clock reads.
            if (counter == 1 ||
                counter * 100 / selected_I.size() != (counter - 1) * 100 / selected_I.size())
                print_progress(iteration_I, counter, selected_I.size());
            // If i generates I, I intersect f^-1(J) = i((f*i)^-1(J)).
            // Work on the free generator basis; no presentation of I is needed.
            auto restricted = measure(profile, &PruningProfile::i_image, [&] {
                return f.generator_lift() * I.generator_map().generator_lift();
            });
            auto inverse = measure(profile, &PruningProfile::i_preimage, [&] {
                if (options.reuse_target_reductions)
                    return pruning_preimage_with_target(restricted, *target_generators,
                        options.projected_preimages);
                return graded_linalg::inverse_image_generators(
                    restricted, f.target()->presentation(),
                    J.generator_map().generator_lift(), options.projected_preimages,
                    true);
            });
            I = measure(profile, &PruningProfile::i_image, [&] {
                return Submodule(I.parent(), I.generator_map().generator_lift() * inverse);
            });
            measure(profile, &PruningProfile::i_reduce, [&] {
                I.minimize_generators(true, true);
            });
        }
        std::cout << std::endl;
        // I decreases, so only the reverse containment needs checking.
        if (measure(profile, &PruningProfile::i_equal, [&] {
                return I.contains(I_);
            }))
            break;
    }
    return {std::move(I), iteration_I};
}

/** Grow K by repeated sums and preimages until a full pass leaves it unchanged. */
template <typename Matrix>
std::pair<graded_linalg::Submodule<Matrix>, int>
compute_pruning_k(const graded_linalg::Homomorphism<Matrix>& eta_X,
                  const std::vector<graded_linalg::Homomorphism<Matrix>>& B_K,
                  PruningProfile* profile,
                  const PruningOptions& options) {
    using Submodule = graded_linalg::Submodule<Matrix>;
    using Hom = graded_linalg::Homomorphism<Matrix>;
    Submodule K = eta_X.domain()->zero_submodule();
    int iteration_K = 0;

    PruningTimer phase(profile ? &profile->k_phase : nullptr);
    for (;;) {
        std::size_t counter = 0;
        iteration_K++;
        if (profile)
            ++profile->k_iterations;
        const Submodule K_ = K;
        auto L = measure(profile, &PruningProfile::k_image, [&] {
            return eta_X.image(K, true, true);
        });
        std::vector<std::size_t> selected_K;
        if (options.adaptive_restriction) {
            selected_K = measure(profile, &PruningProfile::k_adaptive_restriction, [&] {
                auto target = eta_X.target()->presentation();
                target.append_matrix(L.generator_map().generator_lift());
                return pruning_adaptive_restriction(B_K, K.generator_map().generator_lift(), target);
            });
            if (profile) {
                profile->k_restriction_candidates += B_K.size();
                profile->k_restriction_kept += selected_K.size();
            }
        } else {
            selected_K.resize(B_K.size());
            std::iota(selected_K.begin(), selected_K.end(), std::size_t(0));
        }
        int pending = 0;
        // K is fixed throughout this pass. Accumulate image generators
        // and reduce once per bounded batch instead of after every sum.
        auto accumulated = options.amortize_reductions
            ? L.generator_map().generator_lift() : Matrix{};
        auto flush = [&] {
            L = Submodule(L.parent(), std::move(accumulated));
            measure(profile, &PruningProfile::k_reduce, [&] {
                L.minimize_generators(true, true);
            });
            accumulated = L.generator_map().generator_lift();
            pending = 0;
        };
        for (auto map_index : selected_K) {
            const Hom& f = B_K[map_index];
            ++counter;
            // At most one update per percentage, with no clock reads.
            if (counter == 1 ||
                counter * 100 / selected_K.size() != (counter - 1) * 100 / selected_K.size())
                print_progress(iteration_K, counter, selected_K.size());
            if (options.amortize_reductions && K.num_of_gens() == 0) continue;
            auto image = measure(profile, &PruningProfile::k_image, [&] {
                return f.image(K, false, true);
            });
            if (options.amortize_reductions) {
                measure(profile, &PruningProfile::k_sum, [&] {
                    accumulated.append_matrix(image.generator_map().generator_lift());
                });
                // ponytail: at most 16 images between reductions; this
                // bounds batching without introducing a tuning API.
                if (++pending == 16) flush();
            } else {
                L = measure(profile, &PruningProfile::k_sum, [&] {
                    return L.sum(image, true);
                });
            }
        }
        if (options.amortize_reductions && pending) flush();
        // This implementation takes the preimage under eta_X rather than eta_I.
        // Consequently, K is represented as a submodule of X throughout the loop.
        K = measure(profile, &PruningProfile::k_preimage, [&] {
            return eta_X.preimage(L, true, options.projected_preimages,
                                  true);
        });
        std::cout << std::endl;
        // K increases: eta_X(K_) is included in L before taking its preimage.
        if (measure(profile, &PruningProfile::k_equal, [&] {
                return K_.contains(K);
            }))
            break;
    }
    return {std::move(K), iteration_K};
}

/** Coordinate and grid matrices share this high-level module algorithm.
 * profile optionally collects timing and operation counts; nullptr disables it.
 */
template <typename Matrix>
std::pair<graded_linalg::Submodule<Matrix>, graded_linalg::Submodule<Matrix>>
pruning_pair(const graded_linalg::Module<Matrix>& X,
             double epsilon,
             bool quick = false,
             PruningProfile* profile = nullptr,
             const PruningOptions& options = {}) {

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
        auto [eta_X, B_I, composition_k_maps] = prepare_pruning_pair(X, epsilon, profile, options);
        const std::vector<graded_linalg::Homomorphism<Matrix>>& B_K =
            options.composition_enabled() ? composition_k_maps : B_I;

        // Q = 0 imposes no I constraints, but K must still contain ker(eta).
        if (B_I.empty() && B_K.empty()) {
            auto I = X.whole_submodule();
            PruningTimer timer(profile ? &profile->k_preimage : nullptr);
            auto K = eta_X.kernel();
            return {std::move(I), std::move(K)};
        }

        auto [I, iteration_I] = compute_pruning_i(eta_X, B_I, profile, options);
        auto [K, iteration_K] = compute_pruning_k(eta_X, B_K, profile, options);

        std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K"
                  << std::endl;
        std::cout << "Output is " << std::max(iteration_I, iteration_K) * epsilon
                  << "-interleaved with the input" << std::endl;
        return {std::move(I), std::move(K)};
    }
}

template <typename Matrix>
graded_linalg::Module<Matrix> pruning(graded_linalg::Module<Matrix> X,
                                      double epsilon,
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
        if (profile)
            *profile = {};
        PruningTimer total(profile ? &profile->total : nullptr);
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
        result.shift(diagonal_translation(result.presentation(), -epsilon));
        return result;
    }
}

} // namespace stable_decomposition
