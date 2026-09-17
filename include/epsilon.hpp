#pragma once
#include <optional>
#include "types.hpp"

namespace stable_decomposition {
/** Extract epsilon as 1% of the maximum extent of all presentation degrees. */
std::optional<double> calculate_epsilon_from_matrix(const Mat& matrix);
/** Explicit override, otherwise extraction, with fallback epsilon=0.01. */
double get_epsilon(std::optional<double> user_epsilon, const Mat& matrix);
/** Validate epsilon and form the algorithm's diagonal shift (2 epsilon, 2 epsilon). */
graded_linalg::r2degree pruning_shift(double epsilon);
} // namespace stable_decomposition
