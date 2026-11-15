// A small benchmark worker: load outside the timer and send a start event so
// the runner can impose a deadline on the algorithm rather than input loading.
#include "pruning.hpp"
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <thread>

class BenchmarkMemoryGuard {
    std::atomic<bool> finished_{false};
    std::thread thread_;
public:
    explicit BenchmarkMemoryGuard(double limit_mb) {
        if (limit_mb <= 0) return;
        thread_ = std::thread([this, limit_mb] {
            while (!finished_.load()) {
                rusage usage{};
                if (getrusage(RUSAGE_SELF, &usage)) std::_Exit(76);
#ifdef __APPLE__
                const double megabytes = usage.ru_maxrss / (1024.0 * 1024.0);
#else
                const double megabytes = usage.ru_maxrss / 1024.0;
#endif
                if (megabytes > limit_mb) {
                    std::cerr << "BENCHMARK_MEMORY_LIMIT " << megabytes << " MiB\n";
                    std::_Exit(75);
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
            }
        });
    }
    ~BenchmarkMemoryGuard() {
        finished_ = true;
        if (thread_.joinable()) thread_.join();
    }
};

int main(int argc, char** argv) {
    std::ostream report(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    try {
        if (argc != 4 && argc != 5)
            throw std::invalid_argument("Usage: pruning_optimization_worker INPUT EPSILON MASK [RSS_LIMIT_MIB]");
        std::size_t used = 0;
        const std::string epsilon_text(argv[2]), mask_text(argv[3]);
        const double epsilon = std::stod(epsilon_text, &used);
        if (used != epsilon_text.size()) throw std::invalid_argument("Invalid epsilon");
        const auto mask = std::stoul(mask_text, &used);
        if (used != mask_text.size() || mask > 127) throw std::invalid_argument("Invalid optimization mask");
        // Bits follow scaling.SWITCHES; mask 0 is the experiment baseline.
        stable_decomposition::PruningOptions options;
        options.projected_preimages = mask & 1;
        options.amortize_reductions = mask & 2;
        options.full_end_space = mask & 4;
        options.composition_generators = (mask & 8) || (mask & 64);
        options.composition_budget.sparse_nilpotent = mask & 8;
        options.composition_budget.complement_orbit = mask & 64;
        options.adaptive_restriction = mask & 16;
        options.reuse_target_reductions = mask & 32;
        stable_decomposition::Module input(argv[1]);
        input.sort_compatibly();
        BenchmarkMemoryGuard memory(argc == 5 ? std::stod(argv[4]) : 0);
        stable_decomposition::PruningProfile profile;
        const auto start = std::chrono::steady_clock::now();
        report << "START\n" << std::flush;
        const auto result = stable_decomposition::pruning(std::move(input), epsilon, false, &profile, options);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        report << std::setprecision(17) << "RESULT {\"seconds\":" << elapsed
               << ",\"profile_seconds\":" << profile.total.seconds
               << ",\"basis_seconds\":" << profile.basis.seconds
               << ",\"i_iterations\":" << profile.i_iterations
               << ",\"k_iterations\":" << profile.k_iterations
               << ",\"generators_out\":" << result.num_of_gens()
               << ",\"relations_out\":" << result.num_of_rels() << "}\n" << std::flush;
    } catch (const std::exception& error) {
        std::cerr << "Benchmark error: " << error.what() << '\n';
        return 1;
    }
}
