#pragma once
#include "epsilon.hpp"

namespace stable_decomposition {
// Historical spellings. New code should include epsilon.hpp and use epsilon APIs.
std::optional<Epsilon> calculate_delta_from_matrix(const Mat& matrix);
Epsilon get_delta(std::optional<Epsilon> value, const Mat& matrix);
} // namespace stable_decomposition
