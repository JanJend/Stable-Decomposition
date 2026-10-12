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
#include <grlina/hom_interface.hpp>
#include <grlina/mat_csc.hpp>
#include <grlina/mat_grid.hpp>
#include <grlina/mat_coord.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace stable_decomposition {
namespace module_pruning_detail {

// Grid matrices distinguish geometric shifts from stored index degrees.
template <typename Matrix>
auto diagonal_translation(const Matrix&, double amount) ->
    graded_linalg::matrix_geometry_degree_t<Matrix> {
    using Degree = graded_linalg::matrix_geometry_degree_t<Matrix>;
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
    std::is_same_v<typename Matrix::storage_type,
                   graded_linalg::CSCStorage<typename Matrix::index_type>>;

template <class Storage, class Matrix>
auto convert_matrix_storage(const Matrix& matrix) {
    using Index = typename Matrix::index_type;
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
        convert_complex(module.projective_resolution()),
        module.has_complete_projective_resolution()
            ? graded_linalg::ResolutionCompleteness::complete
            : graded_linalg::ResolutionCompleteness::truncated);
    if (module.has_injective_resolution())
        result.set_injective_resolution(convert_complex(module.injective_resolution()));
    return result;
}

} // namespace module_pruning_detail

/** Coordinate and grid matrices share this high-level module algorithm. */
template <typename Matrix>
std::pair<graded_linalg::Submodule<Matrix>, graded_linalg::Submodule<Matrix>>
pruning_pair(const graded_linalg::Module<Matrix>& X,
             double epsilon,
             bool quick = false,
             PruningProfile* profile = nullptr) {

    if constexpr (module_pruning_detail::editable_conversion_required<Matrix>) {
        auto original = X.shared_from_this();
        using EditableStorage = graded_linalg::vec<graded_linalg::vec<typename Matrix::index_type>>;
        auto editable = module_pruning_detail::convert_module_storage<EditableStorage>(X);
        auto parent = std::make_shared<const decltype(editable)>(std::move(editable));
        auto [I, K] =
            stable_decomposition::pruning_pair(*parent, epsilon, quick, profile);
        // Preserve the caller's exact parent identity. Pair results represent
        // embeddings; their optional presentation cache is computed on demand.
        return {graded_linalg::Submodule<Matrix>(
                    original,
                    module_pruning_detail::convert_matrix_storage<typename Matrix::storage_type>(
                        I.generator_map().generator_lift())),
                graded_linalg::Submodule<Matrix>(
                    original,
                    module_pruning_detail::convert_matrix_storage<typename Matrix::storage_type>(
                        K.generator_map().generator_lift()))};
    } else {

        (void)pruning_shift(epsilon); // shared finite/nonnegative validation
        (void)quick; // the high-level algorithm retains the existing ordinary iteration
        using Submodule = graded_linalg::Submodule<Matrix>;
        using Hom = graded_linalg::Homomorphism<Matrix>;
        const auto eps_vec =
            module_pruning_detail::diagonal_translation(X.presentation(), 2 * epsilon);
        Hom eta_X = graded_linalg::Homomorphism<Matrix>::canonical_shift(X, eps_vec);

        auto B = measure(profile, &PruningProfile::basis, [&] {
            return timed_with_progress("End_2eps/0", [&] {
                return graded_linalg::End_2d_0(eta_X, true);
            });
        });

        Submodule I = X.whole_submodule();
        // Preserve the existing empty-family convention: return (X, 0).
        // This differs from initializing K to the shift kernel in the paper.
        if (B.empty())
            return {I, X.zero_submodule()};

        // Repeated intersections shrink I until a full pass leaves it unchanged.
        int iteration_I = 0;
        int counter = 0;
        {
            PruningTimer phase(profile ? &profile->i_phase : nullptr);
            for (;;) {
                counter = 0;
                iteration_I++;
                if (profile)
                    ++profile->i_iterations;
                const Submodule I_ = I;
                auto J = measure(profile, &PruningProfile::i_image, [&] {
                    return eta_X.image(I);
                });
                for (const Hom& f : B) {
                    print_progress(iteration_I, ++counter, B.size());
                    auto preimage = measure(profile, &PruningProfile::i_preimage, [&] {
                        return f.preimage(J);
                    });
                    I = measure(profile, &PruningProfile::i_intersection, [&] {
                        return I.intersection(preimage, false);
                    });
                    measure(profile, &PruningProfile::i_reduce, [&] {
                        I.minimize_generators();
                    });
                }
                print_progress(iteration_I, counter, B.size());
                std::cout << std::endl;
                // I decreases, so only the reverse containment needs checking.
                if (measure(profile, &PruningProfile::i_equal, [&] {
                        return I.contains(I_);
                    }))
                    break;
            }
        }

        // Repeated sums and preimages grow K until a full pass leaves it unchanged.
        Submodule K = X.zero_submodule();
        int iteration_K = 0;

        {
            PruningTimer phase(profile ? &profile->k_phase : nullptr);
            for (;;) {
                counter = 0;
                iteration_K++;
                if (profile)
                    ++profile->k_iterations;
                const Submodule K_ = K;
                auto L = measure(profile, &PruningProfile::k_image, [&] {
                    return eta_X.image(K);
                });
                for (const Hom& f : B) {
                    print_progress(iteration_K, ++counter, B.size());
                    auto image = measure(profile, &PruningProfile::k_image, [&] {
                        return f.image(K);
                    });
                    L = measure(profile, &PruningProfile::k_sum, [&] {
                        return L.sum(image);
                    });
                }
                // This implementation takes the preimage under eta_X rather than eta_I.
                // Consequently, K is represented as a submodule of X throughout the loop.
                K = measure(profile, &PruningProfile::k_preimage, [&] {
                    return eta_X.preimage(L);
                });
                print_progress(iteration_K, counter, B.size());
                std::cout << std::endl;
                // K increases: eta_X(K_) is included in L before taking its preimage.
                if (measure(profile, &PruningProfile::k_equal, [&] {
                        return K_.contains(K);
                    }))
                    break;
            }
        }

        std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K"
                  << std::endl;
        std::cout << "Output is " << std::max(iteration_I, iteration_K) * epsilon
                  << "-interleaved with the input" << std::endl;
        return {I, K};
    }
}

template <typename Matrix>
graded_linalg::Module<Matrix> pruning(graded_linalg::Module<Matrix> X,
                                      double epsilon,
                                      bool quick = false,
                                      PruningProfile* profile = nullptr) {
    if constexpr (module_pruning_detail::editable_conversion_required<Matrix>) {
        using EditableStorage = graded_linalg::vec<graded_linalg::vec<typename Matrix::index_type>>;
        auto editable = module_pruning_detail::convert_module_storage<EditableStorage>(X);
        auto result = stable_decomposition::pruning(
            std::move(editable), epsilon, quick, profile);
        return module_pruning_detail::convert_module_storage<typename Matrix::storage_type>(result);
    } else {
        using Module = graded_linalg::Module<Matrix>;
        if (profile)
            *profile = {};
        PruningTimer total(profile ? &profile->total : nullptr);
        auto parent = std::make_shared<Module>(std::move(X));
        auto [I, K] = stable_decomposition::pruning_pair<Matrix>(*parent, epsilon, quick, profile);
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
        result.shift(module_pruning_detail::diagonal_translation(result.presentation(), -epsilon));
        return result;
    }
}

} // namespace stable_decomposition
