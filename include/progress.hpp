#pragma once
#include <grlina/progress.hpp>

namespace stable_decomposition {
using graded_linalg::timed_with_progress;
void print_progress(int iteration, std::size_t current, std::size_t total);
} // namespace stable_decomposition
