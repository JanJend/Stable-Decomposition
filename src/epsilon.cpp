#include "epsilon.hpp"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "delta.hpp"

namespace stable_decomposition {




std::optional<double> calculate_epsilon_from_matrix(const Mat& M) {
    const auto& col_degrees = M.col_degrees;
    const auto& row_degrees = M.row_degrees;

    if (col_degrees.empty() && row_degrees.empty()) {
        return std::nullopt;
    }

    // Initialize with first available degree
    double min_x, max_x, min_y, max_y;

    if (!col_degrees.empty()) {
        min_x = max_x = col_degrees[0].first;
        min_y = max_y = col_degrees[0].second;
    } else {
        min_x = max_x = row_degrees[0].first;
        min_y = max_y = row_degrees[0].second;
    }

    // Find bounding box of all degrees
    for (const auto& d : col_degrees) {
        min_x = std::min(min_x, d.first);
        max_x = std::max(max_x, d.first);
        min_y = std::min(min_y, d.second);
        max_y = std::max(max_y, d.second);
    }
    for (const auto& d : row_degrees) {
        min_x = std::min(min_x, d.first);
        max_x = std::max(max_x, d.first);
        min_y = std::min(min_y, d.second);
        max_y = std::max(max_y, d.second);
    }

    // Compute max extent
    double extent_x = max_x - min_x;
    double extent_y = max_y - min_y;
    double extent = std::max(extent_x, extent_y);

    if (extent < 1e-10) {
        return std::nullopt;
    }

    return extent * 0.01; // 1% of the extent
}

double get_epsilon(std::optional<double> user_epsilon, const Mat& M) {
    if (user_epsilon.has_value()) {
        (void)pruning_shift(user_epsilon.value());
        std::cout << "Using specified epsilon: " << user_epsilon.value() << std::endl;
        return user_epsilon.value();
    }

    auto calculated = calculate_epsilon_from_matrix(M);
    if (calculated.has_value()) {
        (void)pruning_shift(calculated.value());
        std::cout << "Using calculated epsilon: " << calculated.value() << std::endl;
        return calculated.value();
    }

    std::cout << "Using default epsilon: 0.01" << std::endl;
    return 0.01;
}

// Historical function names retained for source and link compatibility.
std::optional<double> calculate_delta_from_matrix(const Mat& matrix) {
    return calculate_epsilon_from_matrix(matrix);
}
double get_delta(std::optional<double> value, const Mat& matrix) {
    return get_epsilon(value, matrix);
}

} // namespace stable_decomposition
