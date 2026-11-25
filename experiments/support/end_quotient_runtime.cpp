// Compare the explicit two-Hom quotient with the new End_d_0 coefficient route.
#include "epsilon.hpp"
#include "process_memory_guard.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

namespace gl = graded_linalg;
namespace sd = stable_decomposition;
using Clock = std::chrono::steady_clock;
using Coefficients = gl::SparseMatrix<std::ptrdiff_t, gl::DefaultSparseStorage<std::ptrdiff_t>>;

double number(const char* text) {
    std::size_t used = 0;
    const std::string value(text);
    const double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result) || result < 0)
        throw std::invalid_argument("Expected a finite nonnegative number");
    return result;
}

int main(int argc, char** argv) {
    constexpr const char* usage =
        "Usage: end_quotient_benchmark INPUT.scc EPSILON_PERCENT explicit|end-d-0 MEMORY_MB\n"
        "quotient_methods=explicit,end-d-0 hom_algorithm=optimised\n"
        "Internal timings exclude loading, sorting, row-cache construction and shift setup.\n";
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << usage;
        return 0;
    }
    std::ostream report(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    try {
        if (argc != 5) throw std::invalid_argument(usage);
        const double percentage = number(argv[2]), memory_mb = number(argv[4]);
        const std::string method(argv[3]);
        if (method != "explicit" && method != "end-d-0") throw std::invalid_argument(usage);
        sd::ProcessMemoryGuard memory(memory_mb);
        auto module = std::make_shared<sd::Module>(std::string(argv[1]));
        module->sort_compatibly();
        auto source = module->presentation();
        source.compute_rows_forward();
        const auto epsilon = sd::epsilon_from_percent(source, percentage);
        const auto eta = sd::Hom::canonical_shift(module, sd::pruning_shift(epsilon));
        constexpr auto algorithm = gl::HomAlgorithm::optimised;
        report << "START method=" << method << " epsilon=" << sd::format_epsilon(epsilon, 17)
               << " shift=" << sd::format_epsilon(sd::pruning_shift(epsilon), 17) << std::endl;
        const auto start = Clock::now();
        double end0_seconds = 0, seconds = 0;
        std::ptrdiff_t end0_dimension = -1;
        gl::EndQuotientCoefficients<Coefficients> result;
        if (method == "explicit") {
            auto end = gl::hom_space_coefficients(source, source, algorithm);
            end0_seconds = std::chrono::duration<double>(Clock::now() - start).count();
            end0_dimension = end.get_num_cols();
            report << "STAGE End_0_dimension=" << end0_dimension << " seconds=" << end0_seconds << std::endl;
            // This overload performs the old route: shifted Hom, normalization
            // of eta End_0, and ordinary coefficient-space quotient reduction.
            const auto quotient_start = Clock::now();
            result = gl::End_d_0_coefficients(eta, end, algorithm);
            seconds = end0_seconds + std::chrono::duration<double>(Clock::now() - quotient_start).count();
        } else {
            // Deliberately omit E: the supplied-E overload takes the old route.
            result = gl::End_d_0_coefficients(eta, algorithm);
            seconds = std::chrono::duration<double>(Clock::now() - start).count();
        }
        rusage usage_info{};
        if (getrusage(RUSAGE_SELF, &usage_info)) throw std::runtime_error("Cannot read peak RSS");
#ifdef __APPLE__
        const double peak_mb = usage_info.ru_maxrss / (1024.0 * 1024.0);
#else
        const double peak_mb = usage_info.ru_maxrss / 1024.0;
#endif
        report << std::setprecision(17) << "RESULT {\"method\":\"" << method
               << "\",\"hom_algorithm\":\"optimised\",\"seconds\":" << seconds
               << ",\"end0_seconds\":" << end0_seconds
               << ",\"with_E_seconds\":" << seconds - end0_seconds
               << ",\"End_0_dimension\":" << end0_dimension
               << ",\"End_2epsilon_dimension\":" << result.quotient.get_num_cols() + result.subspace.get_num_cols()
               << ",\"quotient_dimension\":" << result.quotient.get_num_cols()
               << ",\"image_dimension\":" << result.subspace.get_num_cols()
               << ",\"peak_rss_mb\":" << peak_mb << "}" << std::endl;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
