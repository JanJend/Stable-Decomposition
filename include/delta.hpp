#pragma once
#include "epsilon.hpp"

namespace stable_decomposition {
// Historical spellings. New code should include epsilon.hpp and use epsilon APIs.
std::optional<double> calculate_delta_from_matrix(const Mat& matrix);
/** Prefer the explicit value, then the heuristic, otherwise 0.01. */
double get_delta(std::optional<double> value, const Mat& matrix);
} // namespace stable_decomposition
