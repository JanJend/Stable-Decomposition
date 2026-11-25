// Build: cmake --build build/local-release --target end_quotient_dimensions
// Run: build/local-release/end_quotient_dimensions INPUT.scc EPSILON_PERCENT
#include "epsilon.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#ifdef __APPLE__
#include <mach/mach.h>
#endif

namespace gl = graded_linalg;
namespace sd = stable_decomposition;
const auto started = std::chrono::steady_clock::now();

void progress(const char* stage) {
    std::cerr << stage << " seconds="
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
#ifdef __APPLE__
    task_vm_info_data_t memory{};
    mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_VM_INFO, reinterpret_cast<task_info_t>(&memory), &count) == KERN_SUCCESS)
        std::cerr << " footprint_GiB=" << memory.phys_footprint / double(1ULL << 30)
                  << " peak_footprint_GiB=" << memory.ledger_phys_footprint_peak / double(1ULL << 30);
#endif
    std::cerr << std::endl;
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("Usage: end_quotient_dimensions INPUT.scc EPSILON_PERCENT");
        std::size_t used = 0;
        const double percent = std::stod(argv[2], &used);
        if (used != std::string(argv[2]).size() || !std::isfinite(percent) || percent < 0)
            throw std::invalid_argument("Epsilon percent must be finite and nonnegative");
        sd::Module module{std::string(argv[1])};
        module.sort_compatibly();
        auto source = module.presentation();
        source.compute_rows_forward();
        const auto shift = sd::pruning_shift(sd::epsilon_from_percent(source, percent));
        std::cerr << std::setprecision(10) << "epsilon_percent=" << percent
                  << " shift=" << sd::format_epsilon(shift) << std::endl;
        auto end = gl::hom_space_coefficients(source, source, gl::HomAlgorithm::optimised);
        std::cout << "End_0_dimension=" << end.get_num_cols() << std::endl;
        progress("End_0 computed");
        std::shared_ptr<const sd::Module> parent = std::make_shared<const sd::Module>(std::move(module));
        const auto eta = sd::Hom::canonical_shift(parent, shift);
        auto quotient = gl::End_d_0_coefficients(eta, end);
        const auto q_dimension = quotient.quotient.get_num_cols();
        std::cout << "End_2e_dimension=" << q_dimension + quotient.subspace.get_num_cols() << std::endl;
        std::cout << "Q_dimension=" << q_dimension << std::endl;
        std::cerr << "coefficient_rows=" << quotient.quotient.get_num_rows()
                  << " rank_eta_End_0=" << quotient.subspace.get_num_cols() << std::endl;
        progress("quotient computed");
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
