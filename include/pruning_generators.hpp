#pragma once

#include <grlina/hom_interface.hpp>
#include <algorithm>
#include <map>
#include <random>
#include <stdexcept>
#include <vector>

namespace stable_decomposition {

// Finite coefficient spaces over F2, rather than graded submodule spans.
// Payloads record coordinates in an already chosen quotient basis.
struct PruningCoefficientBasis {
    using Vector = std::vector<std::size_t>;
    struct Reducer { Vector entries, coordinates; };
    std::map<std::size_t, Reducer> pivots;

    static void add(const Vector& source, Vector& target) {
        graded_linalg::Column_traits<Vector, std::size_t>::add_to(source, target);
    }
    Vector normal(Vector entries, Vector* coordinates = nullptr) const {
        Vector remainder;
        while (!entries.empty()) {
            auto it = pivots.find(entries.back());
            if (it == pivots.end()) {
                remainder.push_back(entries.back());
                entries.pop_back();
            } else {
                add(it->second.entries, entries);
                if (coordinates) add(it->second.coordinates, *coordinates);
            }
        }
        std::reverse(remainder.begin(), remainder.end());
        return remainder;
    }
    bool insert(Vector entries, Vector coordinates = {}) {
        while (!entries.empty()) {
            auto it = pivots.find(entries.back());
            if (it == pivots.end()) {
                const auto pivot = entries.back();
                pivots.emplace(pivot, Reducer{std::move(entries), std::move(coordinates)});
                return true;
            }
            add(it->second.entries, entries);
            add(it->second.coordinates, coordinates);
        }
        return false;
    }
    std::size_t rank() const { return pivots.size(); }
};

template<class Matrix>
PruningCoefficientBasis::Vector pruning_coefficients(const Matrix& matrix) {
    PruningCoefficientBasis::Vector entries;
    for (graded_linalg::matrix_index_t<Matrix> c = 0; c < matrix.get_num_cols(); ++c)
        for (auto r : matrix.column(c))
            entries.push_back(std::size_t(c) * std::size_t(matrix.get_num_rows()) + std::size_t(r));
    return entries;
}

// Z = {F0 -> target relations -> target F0}. Each source column is independent.
// Insert sparse coefficient vectors directly; do not allocate one matrix per lift.
template<class Matrix>
PruningCoefficientBasis pruning_zero_lifts(const Matrix& source, const Matrix& target) {
    using Index = graded_linalg::matrix_index_t<Matrix>;
    const auto source_degrees = graded_linalg::geometric_row_degrees(source);
    const auto relation_degrees = graded_linalg::geometric_col_degrees(target);
    using Degree = typename decltype(source_degrees)::value_type;
    PruningCoefficientBasis zeros;
    for (Index c = 0; c < source.get_num_rows(); ++c)
        for (Index j = 0; j < target.get_num_cols(); ++j)
            if (graded_linalg::Degree_traits<Degree>::smaller_equal(relation_degrees[j], source_degrees[c])) {
                PruningCoefficientBasis::Vector entries;
                for (auto r : target.column(j))
                    entries.push_back(std::size_t(c) * std::size_t(target.get_num_rows()) + std::size_t(r));
                zeros.insert(std::move(entries));
            }
    return zeros;
}

template<class Matrix>
struct PruningMapFamilies {
    std::vector<graded_linalg::Homomorphism<Matrix>> i_maps, k_maps;
    std::size_t original_size = 0, quotient_dimension = 0, end_dimension = 0;
};

/** Q = Hom(M,M(delta)) / eta End(M). Select left End generators for I
 * and right End generators for K. See composition_generators_pruning.tex.
 * The candidates are heuristic; a full-rank certificate makes the result exact.
 * Trusted input: eta is the canonical shift by delta, with identity coefficients.
 */
template<class Matrix>
PruningMapFamilies<Matrix> pruning_composition_generators(
    const graded_linalg::Homomorphism<Matrix>& eta,
    const typename Matrix::degree_type& delta, bool select_composition_generators = true) {
    using Basis = PruningCoefficientBasis;
    using Vector = Basis::Vector;
    PruningMapFamilies<Matrix> result;
    auto source = eta.domain()->presentation();
    if (source.get_num_rows() == 0) return result;
    source.compute_rows_forward();
    const auto& target = eta.target()->presentation();
    // An empty vector returned by a conversion cap is not a zero Hom space.
    // Materialize the actual bases; benchmark process limits bound resources.
    auto hom = graded_linalg::hom_space_basis_new(source, target, true, false, true, 0);
    auto end = graded_linalg::hom_space_basis_new(source, source, true, false, true, 0);
    result.end_dimension = end.size();
    auto quotient = pruning_zero_lifts(source, target);
    // dim L = dim Hom + dim Z. Retain the old empty-lift-family convention,
    // including when Q is zero but the original lift quotient was nonzero.
    const auto shifted_lifts = hom.size() + quotient.rank();
    const auto ordinary_lifts = end.size() + pruning_zero_lifts(source, source).rank();
    if (shifted_lifts < ordinary_lifts)
        throw std::logic_error("Canonical shift decreased the lift-space dimension");
    result.original_size = shifted_lifts - ordinary_lifts;
    if (!result.original_size) return result;
    for (const auto& e : end) quotient.insert(pruning_coefficients(e));

    Basis coordinates;
    std::vector<Matrix> q;
    for (auto& f : hom) {
        auto normal = quotient.normal(pruning_coefficients(f));
        if (coordinates.insert(std::move(normal), {q.size()}))
            q.push_back(std::move(f));
    }
    const auto d = q.size();
    result.quotient_dimension = d;
    if (!d) return result;
    if (!select_composition_generators) {
        for (auto& f : q)
            result.i_maps.emplace_back(eta.domain(), eta.target(), std::move(f));
        return result;
    }

    std::vector<Vector> candidates(1);
    for (std::size_t j = 0; j < d; ++j) candidates.front().push_back(j);
    // ponytail: bounded greedy search, not a minimum-generator solver. Basis
    // candidates guarantee progress; more candidates can improve compactness.
    std::mt19937 random(923407);
    for (int trial = 0; trial < 64; ++trial) {
        Vector candidate;
        for (std::size_t j = 0; j < d; ++j) if (random() & 1) candidate.push_back(j);
        candidates.push_back(std::move(candidate));
    }
    for (std::size_t j = 0; j < d; ++j) candidates.push_back({j});

    for (bool left : {true, false}) {
        // Cache action matrices in Q coordinates, once per side.
        auto acting_end = end;
        if (left) for (auto& e : acting_end) e.shift(delta);
        std::vector<std::vector<Vector>> actions(d);
        for (std::size_t j = 0; j < d; ++j) {
            for (const auto& e : acting_end) {
                auto product = left ? e * q[j] : q[j] * e;
                auto normal = quotient.normal(pruning_coefficients(product));
                Vector value;
                if (!coordinates.normal(std::move(normal), &value).empty())
                    throw std::logic_error("Composition escaped the Hom quotient");
                actions[j].push_back(std::move(value));
            }
        }
        auto add_orbit = [&](const Vector& candidate, Basis& span) {
            for (std::size_t a = 0; a < end.size(); ++a) {
                Vector value;
                for (auto j : candidate) Basis::add(actions[j][a], value);
                span.insert(std::move(value));
            }
        };
        Basis generated;
        auto& selected = left ? result.i_maps : result.k_maps;
        while (generated.rank() < d) {
            std::size_t best_rank = generated.rank();
            const Vector* best = nullptr;
            for (const auto& candidate : candidates) {
                auto trial = generated;
                add_orbit(candidate, trial);
                if (trial.rank() > best_rank) {
                    best_rank = trial.rank();
                    best = &candidate;
                }
                if (best_rank == d) break;
            }
            if (!best) throw std::logic_error("Composition generators failed the rank certificate");
            add_orbit(*best, generated);
            Matrix representative = q[best->front()];
            for (std::size_t k = 1; k < best->size(); ++k)
                q[(*best)[k]].add_matrix_to(representative);
            selected.emplace_back(eta.domain(), eta.target(), std::move(representative));
        }
    }
    return result;
}

/** Select a basis of the restrictions of maps to U, modulo a fixed target
 * submodule C (including the target relations). Each source-generator grade
 * has its own cached scalar reducer: later-born C columns must not be used.
 * This preserves all I constraints, or the K image sum modulo eta(K).
 */
template<class Matrix>
std::vector<std::size_t> pruning_phase_map_basis(
    const std::vector<graded_linalg::Homomorphism<Matrix>>& maps,
    const Matrix& source_generators, const Matrix& target_generators) {
    using Index = graded_linalg::matrix_index_t<Matrix>;
    using Degree = typename Matrix::degree_type;
    using Basis = PruningCoefficientBasis;
    std::vector<std::size_t> selected;
    if (maps.empty() || source_generators.get_num_cols() == 0) return selected;
    const auto source_degrees = graded_linalg::geometric_col_degrees(source_generators);
    const auto target_degrees = graded_linalg::geometric_col_degrees(target_generators);
    std::map<Degree, Basis> reducers;
    for (const auto& degree : source_degrees) {
        auto [it, inserted] = reducers.try_emplace(degree);
        if (!inserted) continue;
        for (Index j = 0; j < target_generators.get_num_cols(); ++j)
            if (graded_linalg::Degree_traits<Degree>::smaller_equal(target_degrees[j], degree)) {
                Basis::Vector entries;
                for (auto r : target_generators.column(j)) entries.push_back(std::size_t(r));
                it->second.insert(std::move(entries));
            }
    }
    Basis span;
    for (std::size_t j = 0; j < maps.size(); ++j) {
        const auto restricted = maps[j].generator_lift() * source_generators;
        Basis::Vector entries;
        for (Index c = 0; c < restricted.get_num_cols(); ++c) {
            Basis::Vector column;
            for (auto r : restricted.column(c)) column.push_back(std::size_t(r));
            auto normal = reducers.at(source_degrees[c]).normal(std::move(column));
            for (auto r : normal)
                entries.push_back(std::size_t(c) * std::size_t(restricted.get_num_rows()) + r);
        }
        if (span.insert(std::move(entries))) selected.push_back(j);
    }
    return selected;
}

// Preserve the exact ambient generator basis: presentation minimization would
// also cancel rows and would therefore be wrong here.
template<class Matrix>
Matrix pruning_reduce_target(Matrix generators, bool dense_pivots) {
    using Index = graded_linalg::matrix_index_t<Matrix>;
    if constexpr (std::tuple_size<typename Matrix::degree_type>::value == 2) {
        auto payload = graded_linalg::empty_matrix_like(generators, generators.get_num_cols(), Index(0));
        payload.col_degrees = generators.col_degrees;
        return graded_linalg::bigraded_reduce(std::move(generators), std::move(payload), true, dense_pivots).image;
    } else {
        auto free_presentation = graded_linalg::empty_matrix_like(generators, Index(0), generators.get_num_rows());
        free_presentation.row_degrees = generators.row_degrees;
        auto parent = std::make_shared<const graded_linalg::Module<Matrix>>(std::move(free_presentation));
        graded_linalg::Submodule<Matrix> image(parent, std::move(generators));
        image.minimize_generators(true, dense_pivots);
        return image.generator_map().generator_lift();
    }
}

template<class Matrix>
Matrix pruning_preimage_with_target(const Matrix& map, const Matrix& target,
                                    bool projected, bool dense_pivots) {
    auto empty = graded_linalg::empty_matrix_like(target, graded_linalg::matrix_index_t<Matrix>(0), target.get_num_rows());
    empty.row_degrees = target.row_degrees;
    return graded_linalg::inverse_image_generators(map, target, empty, projected, dense_pivots);
}

} // namespace stable_decomposition
