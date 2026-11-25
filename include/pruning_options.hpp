#pragma once
#include <grlina/hom_generators.hpp>
#include <algorithm>
#include <thread>

namespace stable_decomposition {

// Optional optimizations default to on, except the experimental updated I target.
// Hom uses Algorithm A-1/2.
// Dense pivots, fused I, deferred K minimization and throttled progress are unconditional.
struct PruningOptions {
    bool projected_preimages = true;
    bool composition_generators = true; // Master switch for both independently enabled methods.
    bool amortize_reductions = true;
    graded_linalg::HomAlgorithm hom_algorithm = graded_linalg::HomAlgorithm::optimised;
    bool adaptive_restriction = true; // Section 1.3: select independent restrictions at each pass start.
    bool reuse_target_reductions = true;
    bool single_solve_quotient = true; // A-1/2: eliminate unshiftable pairs in one Hom solve.
    bool update_i_target = false; // Fast I: refresh eta(I) before each map, rather than once per pass.
    bool minimize_i_target = false; // Redundant shifted generators are valid in the preimage solve.
    bool parallel_i = false; // Fixed-target I only: intersect independent groups of constraints.
    bool parallel_k = false; // Sum independent groups of images before the pass's single preimage.
    graded_linalg::CompositionGeneratorOptions composition_budget;
    bool full_end_space = false; // Benchmark control: omit the quotient by eta End_0.
    // Detect at runtime; fall back to one worker if the core count is unknown.
    int threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    bool composition_enabled() const {
        return composition_generators && (composition_budget.nilpotent_reduction || composition_budget.complement_orbit);
    }
};

} // namespace stable_decomposition
