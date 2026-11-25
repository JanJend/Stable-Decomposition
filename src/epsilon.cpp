#include "epsilon.hpp"
#include "delta.hpp"
#include <algorithm>
#include <iostream>

namespace stable_decomposition {

Epsilon epsilon_from_percent(const Mat& matrix, double percent) {
    if (!std::isfinite(percent) || percent < 0)
        throw std::invalid_argument("epsilon percentage must be finite and nonnegative");
    Epsilon lower{}, upper{};
    bool first = true;
    for (const auto* degrees : {&matrix.row_degrees, &matrix.col_degrees}) {
        for (const auto& degree : *degrees) {
            for (std::size_t axis = 0; axis < 2; ++axis) {
                if (!std::isfinite(degree[axis]))
                    throw std::invalid_argument("Bounding-box coordinates must be finite");
                lower[axis] = first ? degree[axis] : std::min(lower[axis], degree[axis]);
                upper[axis] = first ? degree[axis] : std::max(upper[axis], degree[axis]);
            }
            first = false;
        }
    }
    Epsilon epsilon{};
    for (std::size_t axis = 0; axis < 2; ++axis)
        epsilon[axis] = static_cast<double>((static_cast<long double>(upper[axis]) - lower[axis]) * (percent / 100.0L));
    (void)pruning_shift(epsilon);
    return epsilon;
}

std::optional<Epsilon> calculate_epsilon_from_matrix(const Mat& matrix) {
    if (matrix.row_degrees.empty() && matrix.col_degrees.empty()) return std::nullopt;
    return epsilon_from_percent(matrix, 1);
}

Epsilon get_epsilon(std::optional<Epsilon> user_epsilon, const Mat& matrix, double percent) {
    const auto epsilon = user_epsilon ? *user_epsilon : epsilon_from_percent(matrix, percent);
    (void)pruning_shift(epsilon);
    std::cout << (user_epsilon ? "Using specified epsilon: " : "Using calculated epsilon: ")
              << format_epsilon(epsilon) << std::endl;
    return epsilon;
}

// Historical function names retained as aliases.
std::optional<Epsilon> calculate_delta_from_matrix(const Mat& matrix) {
    return calculate_epsilon_from_matrix(matrix);
}
Epsilon get_delta(std::optional<Epsilon> value, const Mat& matrix) {
    return get_epsilon(value, matrix);
}
} // namespace stable_decomposition
