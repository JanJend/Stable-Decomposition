#pragma once
#include <grlina/hom_generators.hpp>

namespace stable_decomposition {

// Optional optimizations default to on; Hom uses Algorithm A-1/2.
// Dense pivots, fused I, deferred K minimization and throttled progress are unconditional.
struct PruningOptions {
    bool projected_preimages = true;
    bool composition_generators = true; // Master switch for both independently enabled methods.
    bool amortize_reductions = true;
    graded_linalg::HomAlgorithm hom_algorithm = graded_linalg::HomAlgorithm::optimised;
    bool adaptive_restriction = true; // Section 1.3: select independent restrictions at each pass start.
    bool reuse_target_reductions = true;
    graded_linalg::CompositionGeneratorOptions composition_budget;
    bool full_end_space = false; // Benchmark control: omit the quotient by eta End_0.
    bool composition_enabled() const {
        return composition_generators && (composition_budget.sparse_nilpotent || composition_budget.complement_orbit);
    }
};

} // namespace stable_decomposition
