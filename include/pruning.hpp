/** Pruning-specific entry points. Legacy helpers remain available transitively. */
#pragma once
#include "types.hpp"
#include "algebra_compat.hpp"
#include "delta.hpp"
#include "progress.hpp"

namespace stable_decomposition {

/** Compute the pruning pair (I,K); matrix columns use the parent's F0 basis. */
std::pair<Mat, Mat> pruning_pair(Mat& M, double delta, bool quick = false);
std::pair<OwnedSubmodule, OwnedSubmodule> pruning_pair(
    std::shared_ptr<const PModule> module, double delta, bool quick = false);

/** Compute the pruned module I/K. The legacy matrix overload remains available. */
Mat pruning(Mat& M, double delta, bool quick = false);
PModule pruning(PModule module, double delta, bool quick = false);

} // namespace stable_decomposition
