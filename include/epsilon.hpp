#pragma once
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include "types.hpp"

namespace stable_decomposition {
using Epsilon = graded_linalg::r2degree;

/** Validate every component and form the algorithm's shift 2 epsilon. */
template <typename Scalar, std::size_t Dimension>
std::array<Scalar, Dimension> pruning_shift(std::array<Scalar, Dimension> epsilon) {
    for (auto& component : epsilon) {
        const long double doubled = 2.0L * component;
        if (component < 0 || !std::isfinite(doubled) ||
            doubled > std::numeric_limits<Scalar>::max())
            throw std::invalid_argument("epsilon and 2*epsilon must be finite and nonnegative in every coordinate");
        component = static_cast<Scalar>(doubled);
    }
    return epsilon;
}
Epsilon pruning_shift(double epsilon);

template <typename Scalar, std::size_t Dimension>
std::string format_epsilon(const std::array<Scalar, Dimension>& epsilon, int precision = 6) {
    std::ostringstream out;
    out.precision(precision);
    out << '(';
    for (std::size_t i = 0; i < Dimension; ++i) out << (i ? "," : "") << epsilon[i];
    out << ')';
    return out.str();
}

/** Percentage of each bounding-box side across generator and relation degrees.
 * Empty presentations and zero-length sides produce zero components.
 */
Epsilon epsilon_from_percent(const Mat& matrix, double percent);
std::optional<Epsilon> calculate_epsilon_from_matrix(const Mat& matrix);
/** Explicit coordinate vector, otherwise a percentage (default 1%). */
Epsilon get_epsilon(std::optional<Epsilon> user_epsilon, const Mat& matrix, double percent = 1);
} // namespace stable_decomposition
