/** Pruning-specific entry points. Legacy helpers remain available transitively. */
#pragma once
#include "types.hpp"
#include "algebra_compat.hpp"
#include "epsilon.hpp"
#include "progress.hpp"
#include "pruning_options.hpp"

namespace stable_decomposition {

std::pair<Submodule, Submodule> pruning_pair(const Module& module, const Epsilon& epsilon, const bool quick = false,
                                          PruningProfile* profile = nullptr, const PruningOptions& options = {});

Module pruning(Module module, const Epsilon& epsilon, bool quick = false, PruningProfile* profile = nullptr,
               const PruningOptions& options = {});

/** Compute the epsilon-pruning pair using M(2 epsilon) and its structure map.
 * Matrix columns use the parent's F0 basis. Epsilon is finite and nonnegative.
 */
std::pair<Mat, Mat> pruning_pair(Mat& M, const Epsilon& epsilon, bool quick = false);


/** Return (I/K)(-epsilon): the final translation is half the iteration shift.
 * The legacy matrix overload remains available.
 */
Mat pruning(Mat& M, const Epsilon& epsilon, bool quick = false);

// Separate diagnostic copy; the original matrix implementation stays uninstrumented.
Mat pruning_profiled(Mat& M, const Epsilon& epsilon, bool quick, PruningProfile* profile);


// Scalar overloads use the same amount in every coordinate.
std::pair<Submodule, Submodule> pruning_pair(const Module& module, double epsilon, bool quick = false,
                                          PruningProfile* profile = nullptr, const PruningOptions& options = {});
Module pruning(Module module, double epsilon, bool quick = false, PruningProfile* profile = nullptr,
               const PruningOptions& options = {});
std::pair<Mat, Mat> pruning_pair(Mat& M, double epsilon, bool quick = false);
Mat pruning(Mat& M, double epsilon, bool quick = false);
Mat pruning_profiled(Mat& M, double epsilon, bool quick, PruningProfile* profile);

} // namespace stable_decomposition

// Templated module overloads for coordinate and grid matrices.
#include "module_pruning.hpp"
