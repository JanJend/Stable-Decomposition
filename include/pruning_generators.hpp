#pragma once

#include "progress.hpp"
#include <grlina/hom_generators.hpp>
#include <utility>
#include <vector>

namespace stable_decomposition {

/** Find left and right E-module generators of Q = Hom_0(X,X(delta))/(eta E),
 * where E = End_0(X) and eta: X -> X(delta) is the canonical shift map.
 * End_d_0 first supplies representatives f_j of a vector-space basis of Q.
 * Composition gives the actions e * [f] = [e(delta) o f] and [f] * e = [f o e].
 * Sparse-nilpotent compression uses lowering actions and small sparse seeds;
 * complement-orbit compression then tries sums of the remaining nonpivot maps. A nilpotent-ideal certificate permits stopping at any point;
 * cached quotient coordinates replace the full pairwise action table.
 * Return (left generators for I, right generators for K) as Homomorphism
 * representatives; their classes generate Q. profile only records dimensions
 * and family sizes and does not affect the algebra. With full_end_space,
 * use H itself and omit W; composition selection then generates H.
 */
template<class Matrix>
std::pair<std::vector<graded_linalg::Homomorphism<Matrix>>,
          std::vector<graded_linalg::Homomorphism<Matrix>>>
pruning_composition_generators(
    const graded_linalg::Homomorphism<Matrix>& eta,
    const typename Matrix::degree_type& delta,
    graded_linalg::HomAlgorithm hom_algorithm = graded_linalg::HomAlgorithm::optimised,
    PruningProfile* profile = nullptr, bool full_end_space = false,
    const graded_linalg::CompositionGeneratorOptions& generator_options = {}) {
    std::vector<graded_linalg::Homomorphism<Matrix>> endomorphism_basis =
        graded_linalg::module_endomorphism_basis(eta.domain(), hom_algorithm);
    std::vector<graded_linalg::Homomorphism<Matrix>> quotient_basis = full_end_space
        ? graded_linalg::module_hom_space_basis(eta.domain(), eta.target(), hom_algorithm)
        : graded_linalg::End_d_0(eta, endomorphism_basis, false, hom_algorithm);

    // W = span{eta o e_i}. End_d_0 has already chosen Q = H/W; pass W again
    // so the action routine tests equality in that same Q when finding coordinates.
    std::vector<graded_linalg::Homomorphism<Matrix>> quotient_subspace;
    if (!full_end_space)
        for (const graded_linalg::Homomorphism<Matrix>& endomorphism : endomorphism_basis)
            quotient_subspace.push_back(endomorphism.compose(eta));

    auto generators = graded_linalg::budgeted_composition_generators(
        quotient_basis, endomorphism_basis, quotient_subspace, generator_options);
    if (profile) {
        profile->end_dimension = endomorphism_basis.size();
        profile->original_maps = quotient_basis.size();
        profile->quotient_dimension = full_end_space ? 0 : quotient_basis.size();
        profile->i_maps = generators.first.size();
        profile->k_maps = generators.second.size();
    }
    return generators;
}

/** Select independent restrictions to the current generator basis, modulo
 * the fixed target submodule. Grade-aware normalization and independence
 * testing are shared algebra operations in Persistence-Algebra.
 */
template<class Matrix>
std::vector<std::size_t> pruning_adaptive_restriction(
    const std::vector<graded_linalg::Homomorphism<Matrix>>& maps,
    const Matrix& source_generators, const Matrix& target_generators) {
    if (maps.empty() || source_generators.get_num_cols() == 0) return {};
    // These matrices are lifts on the free generator basis of the submodule.
    std::vector<Matrix> restriction_lifts;
    for (const graded_linalg::Homomorphism<Matrix>& map : maps)
        restriction_lifts.push_back(map.generator_lift() * source_generators);
    return graded_linalg::independent_matrix_family_indices_modulo_target_relations(
        std::move(restriction_lifts), target_generators);
}

// Preserve the exact ambient generator basis: presentation minimization would
// also cancel rows and would therefore be wrong here.
template<class Matrix>
Matrix pruning_reduce_target(Matrix generators) {
    if constexpr (std::tuple_size<typename Matrix::degree_type>::value == 2) {
        Matrix payload = graded_linalg::empty_matrix_like(generators, generators.get_num_cols(), 0);
        payload.col_degrees = generators.col_degrees;
        return graded_linalg::bigraded_reduce(std::move(generators), std::move(payload), true, true).image;
    } else {
        Matrix free_presentation = graded_linalg::empty_matrix_like(generators, 0, generators.get_num_rows());
        free_presentation.row_degrees = generators.row_degrees;
        std::shared_ptr<const graded_linalg::Module<Matrix>> parent =
            std::make_shared<const graded_linalg::Module<Matrix>>(std::move(free_presentation));
        graded_linalg::Submodule<Matrix> image(parent, std::move(generators));
        image.minimize_generators(true, true);
        return image.generator_map().generator_lift();
    }
}

template<class Matrix>
Matrix pruning_preimage_with_target(const Matrix& map, const Matrix& target,
                                    bool projected) {
    Matrix empty = graded_linalg::empty_matrix_like(target, 0, target.get_num_rows());
    empty.row_degrees = target.row_degrees;
    return graded_linalg::inverse_image_generators(map, target, empty, projected, true);
}

} // namespace stable_decomposition
