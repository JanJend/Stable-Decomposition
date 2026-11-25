#pragma once

#include "progress.hpp"
#include <grlina/hom_generators.hpp>
#include <utility>
#include <vector>

namespace stable_decomposition {

using PruningCoefficients = graded_linalg::SparseMatrix<std::ptrdiff_t,
    graded_linalg::DefaultSparseStorage<std::ptrdiff_t>>;

// Store common modules once, rather than grades and empty columns for every map.
// A returned map owns its temporary lift; concurrent readers do not share caches.
template<class Matrix>
class PruningMapFamily {
    PruningCoefficients coefficients_;
    std::shared_ptr<const graded_linalg::Module<Matrix>> domain_, target_;
public:
    PruningMapFamily(PruningCoefficients coefficients,
                     std::shared_ptr<const graded_linalg::Module<Matrix>> domain,
                     std::shared_ptr<const graded_linalg::Module<Matrix>> target)
        : coefficients_(std::move(coefficients)), domain_(std::move(domain)), target_(std::move(target)) {
        if (!domain_ || !target_) throw std::invalid_argument("Map family requires common modules");
        const auto rows = target_->num_of_gens(), columns = domain_->num_of_gens();
        if (rows && std::size_t(columns) > std::size_t(std::numeric_limits<std::ptrdiff_t>::max() / rows))
            throw std::overflow_error("Flattened map coordinates exceed ptrdiff_t");
        if (!empty() && coefficients_.get_num_rows() != std::ptrdiff_t(rows) * columns)
            throw std::invalid_argument("Map family has incompatible coefficient coordinates");
    }
    std::size_t size() const { return coefficients_.get_num_cols(); }
    bool empty() const { return size() == 0; }
    graded_linalg::Homomorphism<Matrix> operator[](std::size_t index) const {
        auto lift = graded_linalg::empty_matrix_like(target_->presentation(), domain_->num_of_gens(), target_->num_of_gens());
        graded_linalg::set_geometric_degrees(lift,
            graded_linalg::geometric_row_degrees(domain_->presentation()),
            graded_linalg::geometric_row_degrees(target_->presentation()));
        const auto rows = target_->num_of_gens();
        for (auto entry : coefficients_.column(index)) lift.append_entry(entry / rows, entry % rows);
        return {domain_, target_, std::move(lift)};
    }
};

/** Find left and right E-module generators of Q = Hom_0(X,X(delta))/(eta E),
 * where E = End_0(X) and eta: X -> X(delta) is the canonical shift map.
 * End_d_0 first supplies representatives f_j of a vector-space basis of Q.
 * Composition gives the actions e * [f] = [e(delta) o f] and [f] * e = [f o e].
 * Nilpotent reduction first certifies elementary lowering factorizations,
 * then samples lowering actions. Small ordinary seeds are skipped when direct
 * certificates fill the nilpotent allowance. Complement-orbit then tries sums
 * of remaining nonpivot maps. A nilpotent-ideal certificate permits stopping
 * at any point; shared quotient coordinates are prepared only when needed.
 * Return coefficient columns for the left I and right K generating families;
 * their classes generate Q. profile only records dimensions
 * and family sizes and does not affect the algebra. With full_end_space,
 * use H itself and omit W; composition selection then generates H.
 */
template<class Matrix>
std::pair<PruningCoefficients, PruningCoefficients>
pruning_composition_coefficients(
    const graded_linalg::Homomorphism<Matrix>& eta,
    const typename Matrix::degree_type& delta,
    graded_linalg::HomAlgorithm hom_algorithm = graded_linalg::HomAlgorithm::optimised,
    PruningProfile* profile = nullptr, bool full_end_space = false,
    const graded_linalg::CompositionGeneratorOptions& generator_options = {},
    bool single_solve_quotient = true) {
    auto source = eta.domain()->presentation();
    source.compute_rows_forward();
    const auto& target = eta.target()->presentation();
    auto endomorphism_basis = measure(profile, &PruningProfile::end0, [&] {
        return graded_linalg::hom_space_coefficients(source, source, hom_algorithm);
    });
    using Coefficients = decltype(endomorphism_basis);
    // Keep E for the actions; A-1/2 computes Q with the single-solve routine.
    // Algorithms A/B still use two Hom solves, so reuse E for those.
    auto quotient = measure(profile, &PruningProfile::hom_quotient, [&] { return full_end_space
        ? graded_linalg::EndQuotientCoefficients<Coefficients>{
            graded_linalg::hom_space_coefficients(source, target, hom_algorithm),
            Coefficients(0, endomorphism_basis.get_num_rows())}
        : single_solve_quotient && hom_algorithm == graded_linalg::HomAlgorithm::optimised
            ? graded_linalg::End_d_0_coefficients(eta, hom_algorithm)
            : graded_linalg::End_d_0_coefficients(eta, endomorphism_basis, hom_algorithm); });
    auto selected = measure(profile, &PruningProfile::composition, [&] {
        return graded_linalg::budgeted_composition_generators(
            quotient.quotient, endomorphism_basis, quotient.subspace, source, target, generator_options);
    });
    if (profile) {
        profile->end_dimension = endomorphism_basis.get_num_cols();
        profile->original_maps = quotient.quotient.get_num_cols();
        profile->quotient_dimension = full_end_space ? 0 : quotient.quotient.get_num_cols();
        profile->i_maps = selected.first.get_num_cols();
        profile->k_maps = selected.second.get_num_cols();
    }
    return selected;
}

// Compatibility helper for callers that explicitly need materialized maps.
template<class Matrix>
std::pair<std::vector<graded_linalg::Homomorphism<Matrix>>,
          std::vector<graded_linalg::Homomorphism<Matrix>>>
pruning_composition_generators(
    const graded_linalg::Homomorphism<Matrix>& eta,
    const typename Matrix::degree_type& delta,
    graded_linalg::HomAlgorithm hom_algorithm = graded_linalg::HomAlgorithm::optimised,
    PruningProfile* profile = nullptr, bool full_end_space = false,
    const graded_linalg::CompositionGeneratorOptions& generator_options = {},
    bool single_solve_quotient = true) {
    auto selected = pruning_composition_coefficients(eta, delta, hom_algorithm, profile,
        full_end_space, generator_options, single_solve_quotient);
    return measure(profile, &PruningProfile::map_lifts, [&] { return std::make_pair(
        graded_linalg::homomorphisms_from_coefficients(selected.first, eta.domain(), eta.target()),
        graded_linalg::homomorphisms_from_coefficients(selected.second, eta.domain(), eta.target())); });
}

/** Select independent restrictions to the current generator basis, modulo
 * the fixed target submodule. Grade-aware normalization and independence
 * testing are shared algebra operations in Persistence-Algebra.
 */
template<class Matrix, class MapFamily>
std::vector<std::size_t> pruning_adaptive_restriction(
    const MapFamily& maps,
    const Matrix& source_generators, const Matrix& target_generators, int threads = 1) {
    if (maps.empty() || source_generators.get_num_cols() == 0) return {};
    const auto reducers = graded_linalg::prepare_target_relation_reducers(
        target_generators, graded_linalg::geometric_col_degrees(source_generators), threads);
    graded_linalg::CompositionCoefficientBasis span;
    std::vector<std::size_t> selected;
    // ponytail: 32 temporary restrictions per phase; retain only sparse pivots.
    // Reduce in original order across batches to preserve the selected indices.
    constexpr std::size_t batch_size = 32;
    for (std::size_t begin = 0; begin < maps.size(); begin += batch_size) {
        const auto count = std::min(batch_size, maps.size() - begin);
        std::vector<graded_linalg::CompositionCoefficientBasis::Vector> columns(count);
        graded_linalg::parallel_for(count, threads, [&](std::size_t i) {
            auto restriction = maps[begin + i].generator_lift() * source_generators;
            graded_linalg::normalize_map_modulo_target_relations(restriction, reducers);
            for (std::ptrdiff_t c = 0; c < restriction.get_num_cols(); ++c)
                for (auto row : restriction.column(c))
                    columns[i].push_back(std::size_t(c) * restriction.get_num_rows() + row);
        });
        for (std::size_t i = 0; i < count; ++i)
            if (span.insert(std::move(columns[i]))) selected.push_back(begin + i);
    }
    return selected;
}

// In one fixed parent, equality of stored generators is a sufficient certificate.
// Different representatives conservatively invalidate the cached target.
template<class Matrix>
bool pruning_same_generators(const Matrix& a, const Matrix& b) {
    if (a.get_num_rows() != b.get_num_rows() || a.get_num_cols() != b.get_num_cols() ||
        graded_linalg::geometric_col_degrees(a) != graded_linalg::geometric_col_degrees(b)) return false;
    for (std::ptrdiff_t c = 0; c < a.get_num_cols(); ++c) {
        const auto& left = a.column(c);
        const auto& right = b.column(c);
        if (!std::equal(left.begin(), left.end(), right.begin(), right.end())) return false;
    }
    return true;
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
