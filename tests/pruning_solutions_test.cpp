#include "pruning.hpp"
#include <grlina/isomorphism_test.hpp>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

int main() {
    namespace sd = stable_decomposition;
    const std::filesystem::path directory(PRUNING_SOLUTIONS_DIR);
    // Historical _pru1 fixtures use iteration shift 1, hence epsilon = 0.5.
    constexpr double epsilon = 0.5;
    int failures = 0;
    for (int i = 1; i <= 4; ++i) {
        const auto name = "test" + std::to_string(i);
        try {
            const sd::Module input((directory / (name + ".scc")).string());
            const sd::Module expected((directory / (name + "_pru1.scc")).string());
            const auto result = sd::pruning(input, epsilon);
            // Minimize both presentations inside the exact isomorphism test.
            const bool matches = graded_linalg::is_isomorphic(
                result.presentation(), expected.presentation());
            std::cout << (matches ? "PASS: " : "FAIL: ") << name
                      << " (epsilon=" << epsilon << ", isomorphism)\n";
            if (!matches) ++failures;
        } catch (const std::exception& error) {
            std::cerr << "ERROR: " << name << ": " << error.what() << '\n';
            ++failures;
        }
    }
    return failures ? 1 : 0;
}
