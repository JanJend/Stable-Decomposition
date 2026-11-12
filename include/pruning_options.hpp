#pragma once

namespace stable_decomposition {

// All false preserves the original module algorithm for timing comparisons.
struct PruningOptions {
    bool fuse_i_preimage = false;
    bool defer_k_image_minimization = false;
    bool projected_preimages = false;
    bool dense_kernel_pivots = false;
    bool throttle_progress = false;
    bool composition_generators = false;
    bool amortize_reductions = false;
    bool hom_classes = false;
    bool phase_quotient = false;
    bool reuse_target_reductions = false;
};

} // namespace stable_decomposition
