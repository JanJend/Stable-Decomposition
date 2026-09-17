/** Pruning-specific entry points. Legacy helpers remain available transitively. */
#pragma once
#include "types.hpp"
#include "algebra_compat.hpp"
#include "epsilon.hpp"
#include "progress.hpp"

namespace stable_decomposition {

using Module = graded_linalg::R2Module<index_t>;
using Submodule = graded_linalg::Submodule<Mat>;

/** Module-based pruning pair (I, K), with I a submodule of module and K a
 * submodule of I. Framework only; throws std::logic_error until implemented.
 */
std::pair<Submodule, Submodule> pruning_pair(const Module& module, double epsilon);

/** Compute the epsilon-pruning pair using M(2 epsilon) and its structure map.
 * Matrix columns use the parent's F0 basis. Epsilon is finite and nonnegative.
 */
std::pair<Mat, Mat> pruning_pair(Mat& M, double epsilon, bool quick = false);
std::pair<Submodule, Submodule> pruning_pair(
    std::shared_ptr<const Module> module, double epsilon, bool quick = false);

/** Return (I/K)(-epsilon): the final translation is half the iteration shift.
 * The legacy matrix overload remains available.
 */
Mat pruning(Mat& M, double epsilon, bool quick = false);
Module pruning(Module module, double epsilon, bool quick = false);

} // namespace stable_decomposition
