#pragma once
#include <optional>
#include "types.hpp"

namespace stable_decomposition {
/** One percent of the maximum coordinate extent over all presentation degrees.
 * Empty or effectively point-supported degree sets yield nullopt.
 */
std::optional<double> calculate_delta_from_matrix(const Mat& matrix);
/** Prefer the explicit value, then the heuristic, otherwise 0.01. */
double get_delta(std::optional<double> user_delta, const Mat& matrix);
} // namespace stable_decomposition
