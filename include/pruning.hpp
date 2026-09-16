/** Pruning-specific entry points. Legacy helpers remain available transitively. */
#pragma once
#include "types.hpp"
#include "algebra_compat.hpp"
#include "epsilon.hpp"
#include "progress.hpp"

namespace stable_decomposition {

/** Compute the epsilon-pruning pair using M(2 epsilon) and its structure map.
 * Matrix columns use the parent's F0 basis. Epsilon is finite and nonnegative.
 */
std::pair<Mat, Mat> pruning_pair(Mat& M, double epsilon, bool quick = false);
std::pair<OwnedSubmodule, OwnedSubmodule> pruning_pair(
    std::shared_ptr<const PModule> module, double epsilon, bool quick = false);

/** Return (I/K)(-epsilon): the final translation is half the iteration shift.
 * The legacy matrix overload remains available.
 */
Mat pruning(Mat& M, double epsilon, bool quick = false);
PModule pruning(PModule module, double epsilon, bool quick = false);

} // namespace stable_decomposition
