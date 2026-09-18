// Compare matrix and module pruning, including an exact isomorphism decision.
#include "pruning.hpp"
#include <grlina/isomorphism_test.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace sd = stable_decomposition;
using DegreeCounts = std::map<std::pair<double, double>, std::size_t>;

namespace {

DegreeCounts degree_counts(const graded_linalg::vec<graded_linalg::r2degree>& degrees) {
    DegreeCounts counts;
    for (const auto& d : degrees) ++counts[{d.first, d.second}];
    return counts;
}

bool compare_degrees(const char* label, const DegreeCounts& old_counts,
                     const DegreeCounts& new_counts) {
    const bool equal = old_counts == new_counts;
    std::cout << label << ": " << (equal ? "PASS" : "FAIL") << '\n';
    if (equal) return true;
    auto all = old_counts;
    all.insert(new_counts.begin(), new_counts.end());
    std::size_t shown = 0;
    for (const auto& entry : all) {
        auto old_it = old_counts.find(entry.first), new_it = new_counts.find(entry.first);
        const auto old_n = old_it == old_counts.end() ? 0 : old_it->second;
        const auto new_n = new_it == new_counts.end() ? 0 : new_it->second;
        if (old_n != new_n && shown++ < 8)
            std::cout << "  degree (" << entry.first.first << ", " << entry.first.second
                      << "): old multiplicity=" << old_n << ", new=" << new_n << '\n';
    }
    if (shown > 8) std::cout << "  ... " << shown - 8 << " more differing degrees\n";
    return false;
}

struct Axis {
    std::vector<double> values;
    bool sampled = false;
};

Axis make_axis(const sd::Mat& old_p, const sd::Mat& new_p, bool x_axis) {
    std::vector<double> values;
    for (const auto* p : {&old_p, &new_p}) {
        for (const auto* degrees : {&p->row_degrees, &p->col_degrees}) {
            for (const auto& d : *degrees) {
                const double value = x_axis ? d.first : d.second;
                if (!std::isfinite(value))
                    throw std::invalid_argument("Comparison grid requires finite presentation degrees");
                values.push_back(value);
            }
        }
    }
    if (values.empty()) values.push_back(0.0);
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    // Include the region before all births and a point beyond all changes.
    const double before = std::nextafter(values.front(), -std::numeric_limits<double>::infinity());
    const double after = std::nextafter(values.back(), std::numeric_limits<double>::infinity());
    if (std::isfinite(before)) values.insert(values.begin(), before);
    if (std::isfinite(after)) values.push_back(after);

    // Keep large inputs usable; explicitly report when the grid is sampled.
    constexpr std::size_t limit = 64;
    if (values.size() <= limit) return {std::move(values), false};
    std::vector<double> sampled;
    for (std::size_t i = 0; i < limit; ++i)
        sampled.push_back(values[i * (values.size() - 1) / (limit - 1)]);
    return {std::move(sampled), true};
}

bool compare_file(const std::filesystem::path& path, double epsilon, bool quick) {
    std::cout << "\n=== " << path << " | epsilon=" << epsilon
              << " | quick=" << std::boolalpha << quick << " ===\n" << std::flush;
    sd::Module loaded(path.string());
    sd::Mat input = loaded.presentation();
    input.sort_compatibly();
    input.validate();
    std::cout << "Input: " << input.get_num_rows() << " generators, "
              << input.get_num_cols() << " relations\n";

    // The matrix overload mutates its argument. Both paths start from independent
    // copies of exactly the same sorted presentation, without higher resolutions.
    sd::Mat old_input = input;
    std::cout << "Running matrix pruning...\n" << std::flush;
    sd::PruningProfile old_profile, new_profile;
    auto start = std::chrono::steady_clock::now();
    sd::Module old_result(sd::pruning_profiled(old_input, epsilon, quick, &old_profile));
    const double old_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "\nRunning module pruning...\n" << std::flush;
    start = std::chrono::steady_clock::now();
    sd::Module new_result = sd::pruning(sd::Module(input), epsilon, quick, &new_profile);
    const double new_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "\nRuntime: old=" << old_seconds << " s, new=" << new_seconds << " s\n";
    sd::print_pruning_comparison(std::cout, old_profile, new_profile);
    old_result.presentation().validate();
    new_result.presentation().validate();
    std::cout << "Presentation validation: PASS\n";

    const auto& old_p = old_result.presentation();
    const auto& new_p = new_result.presentation();
    const bool identical = old_p.row_degrees == new_p.row_degrees &&
                           old_p.col_degrees == new_p.col_degrees && old_p.data == new_p.data;
    std::cout << "Exact output matrices: " << (identical ? "identical" : "different (diagnostic only)") << '\n';

    // Normalize copies with the same minimizer before comparing graded Betti
    // degrees. This avoids treating redundant presentation data as an invariant.
    sd::Module old_min = old_result, new_min = new_result;
    old_min.minimize();
    new_min.minimize();
    std::cout << "Minimized generators: old=" << old_min.number_of_generators()
              << ", new=" << new_min.number_of_generators()
              << "; relations: old=" << old_min.number_of_relations()
              << ", new=" << new_min.number_of_relations() << '\n';
    bool pass = compare_degrees("Generator degree multiset (beta_0)",
        degree_counts(old_min.presentation().row_degrees), degree_counts(new_min.presentation().row_degrees));
    const bool relations_equal = compare_degrees("Relation degree multiset (beta_1)",
        degree_counts(old_min.presentation().col_degrees), degree_counts(new_min.presentation().col_degrees));
    pass = pass && relations_equal;

    // Evaluate the original outputs directly, independently of normalization.
    const auto xs = make_axis(old_p, new_p, true), ys = make_axis(old_p, new_p, false);
    std::size_t mismatches = 0;
    for (double x : xs.values) {
        for (double y : ys.values) {
            const auto old_dim = old_result.dimension_at({x, y});
            const auto new_dim = new_result.dimension_at({x, y});
            if (old_dim != new_dim && mismatches++ < 8)
                std::cout << "  Fibre mismatch at (" << x << ", " << y
                          << "): old=" << old_dim << ", new=" << new_dim << '\n';
        }
    }
    std::cout << "Hilbert function: " << (mismatches == 0 ? "PASS" : "FAIL")
              << " (" << mismatches << " mismatches / " << xs.values.size() * ys.values.size()
              << " points; " << (xs.sampled || ys.sampled ? "sampled" : "full common critical") << " grid)\n";
    pass = pass && mismatches == 0;
    start = std::chrono::steady_clock::now();
    const bool isomorphic = graded_linalg::is_isomorphic(
        old_min.presentation(), new_min.presentation(), true);
    const double iso_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "Isomorphism: " << (isomorphic ? "PASS" : "FAIL") << " (" << iso_seconds << " s)\n";
    pass = pass && isomorphic;
    std::cout << (pass ? "PASS: pruning outputs are isomorphic and all invariants agree.\n"
                      : "FAIL: pruning output comparison failed.\n");
    return pass;
}

} // namespace

int main(int argc, char** argv) {
    try {
        double epsilon = 0.5;
        bool quick = false;
        std::vector<std::filesystem::path> paths;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                std::cout << "Usage: compare_pruning [--epsilon VALUE] [--quick] [FILE ...]\n"
                          << "Defaults: epsilon=0.5, quick=false; toy_example_1.scc, toy_example_2.scc,\n"
                          << "and no_columns_test.scc from " << PRUNING_TEST_PRESENTATIONS_DIR << "\n"
                          << "Bare filenames are also looked up in that directory.\n"
                          << "Exit status: 0=isomorphic and invariants agree, 1=mismatch, 2=invalid arguments or runtime error.\n"
                          << "--quick is passed to both overloads; the module implementation currently ignores it.\n";
                return 0;
            } else if (arg == "--epsilon") {
                if (++i == argc) throw std::invalid_argument("--epsilon requires a value");
                std::size_t consumed = 0;
                epsilon = std::stod(argv[i], &consumed);
                if (consumed != std::string(argv[i]).size()) throw std::invalid_argument("Invalid epsilon");
            } else if (arg == "--quick") {
                quick = true;
            } else if (!arg.empty() && arg.front() == '-') {
                throw std::invalid_argument("Unknown option: " + arg);
            } else {
                paths.emplace_back(arg);
            }
        }
        if (!std::isfinite(epsilon) || epsilon < 0 || !std::isfinite(2 * epsilon))
            throw std::invalid_argument("epsilon and 2*epsilon must be finite and nonnegative");
        if (paths.empty()) paths = {"toy_example_1.scc", "toy_example_2.scc", "no_columns_test.scc"};
        std::cout << std::setprecision(17);
        std::size_t passed = 0, failed = 0, errors = 0;
        for (auto path : paths) {
            if (!std::filesystem::exists(path) && !path.has_parent_path())
                path = std::filesystem::path(PRUNING_TEST_PRESENTATIONS_DIR) / path;
            try {
                if (compare_file(path, epsilon, quick)) ++passed;
                else ++failed;
            } catch (const std::exception& error) {
                ++errors;
                std::cerr << "ERROR for " << path << ": " << error.what() << '\n';
            }
        }
        std::cout << "\nSummary: " << passed << " passed, " << failed << " mismatched, " << errors << " errors.\n";
        return errors ? 2 : failed ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 2;
    }
}
