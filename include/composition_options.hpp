#pragma once
#include <grlina/hom_generators.hpp>
#include <ostream>

namespace stable_decomposition {

// Shared by the pruning CLI, the Hom benchmark parser, and profile output.
struct CompositionCountOption {
    const char* name;
    std::size_t graded_linalg::CompositionGeneratorOptions::*member;
    const char* help;
};
inline constexpr CompositionCountOption composition_count_options[] = {
    {"composition-orbit-rounds", &graded_linalg::CompositionGeneratorOptions::orbit_rounds,
     "Maximum orbit attempts per side (default: 12; 0 disables)"},
    {"composition-orbit-min-gain", &graded_linalg::CompositionGeneratorOptions::orbit_min_gain,
     "Net maps saved for a productive attempt (default: 1; zero savings always count as low)"},
    {"composition-orbit-patience", &graded_linalg::CompositionGeneratorOptions::orbit_patience,
     "Stop after this many consecutive low-gain attempts (default: 3; 0 disables early stopping)"},
    {"composition-orbit-random-seed", &graded_linalg::CompositionGeneratorOptions::orbit_random_seed,
     "Reproducible orbit seed (default: 923407)"},
    {"composition-orbit-pair-visits", &graded_linalg::CompositionGeneratorOptions::orbit_max_pair_visits,
     "Total supported orbit pairs per side (default: 0 = rounds * dim End)"},
    {"composition-orbit-seed-terms", &graded_linalg::CompositionGeneratorOptions::orbit_max_seed_terms,
     "Input entries per orbit seed (default: 65536)"},
    {"composition-orbit-product-terms", &graded_linalg::CompositionGeneratorOptions::orbit_max_product_terms,
     "Symbolic terms per orbit product (default: 65536)"},
    {"composition-nilpotent-rounds", &graded_linalg::CompositionGeneratorOptions::nilpotent_rounds,
     "Maximum nilpotent sampling rounds (default: 8)"},
    {"composition-sparse-seeds", &graded_linalg::CompositionGeneratorOptions::max_seeds,
     "Maximum sparse seed attempts per side (default: 8)"},
    {"composition-sparse-products", &graded_linalg::CompositionGeneratorOptions::max_products,
     "Sparse-method product budget per side (default: 8192)"},
    {"composition-sparse-pair-visits", &graded_linalg::CompositionGeneratorOptions::max_pair_visits,
     "Sparse-method pair-visit budget per side (default: 131072)"},
    {"composition-sparse-seed-basis", &graded_linalg::CompositionGeneratorOptions::max_seed_basis,
     "Basis summands per sparse seed (default: 8)"},
    {"composition-sparse-seed-terms", &graded_linalg::CompositionGeneratorOptions::max_seed_terms,
     "Input entries per sparse seed (default: 256)"},
    {"composition-sparse-product-terms", &graded_linalg::CompositionGeneratorOptions::max_product_terms,
     "Symbolic terms per sparse product (default: 65536)"},
};

inline void print_composition_options(std::ostream& out, const graded_linalg::CompositionGeneratorOptions& options) {
    out << "Composition generator search: composition-orbit-seed-mode=" << options.orbit_seed_mode;
    for (const auto& option : composition_count_options)
        out << ' ' << option.name << '=' << options.*(option.member);
    out << '\n';
}

// The standalone benchmark has a small hand-written parser. Reject signs,
// fractions and overflow here rather than silently converting to an unsigned count.
inline bool parse_composition_count(const std::string& name, const std::string& value,
                                    graded_linalg::CompositionGeneratorOptions& options) {
    for (const auto& option : composition_count_options) if (name == "--" + std::string(option.name)) {
        if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
            throw std::invalid_argument(name + " requires a nonnegative integer");
        const auto count = std::stoull(value);
        if (count > std::numeric_limits<std::size_t>::max()) throw std::out_of_range(name);
        options.*(option.member) = static_cast<std::size_t>(count);
        return true;
    }
    return false;
}

} // namespace stable_decomposition
