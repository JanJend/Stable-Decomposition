// Trace the production pruning path; stdout is flushed after every timer event.
#include "pruning.hpp"
#include "composition_options.hpp"
#include "process_memory_guard.hpp"
#include <cmath>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <map>
#include <sys/resource.h>

namespace sd = stable_decomposition;
using Clock = std::chrono::steady_clock;
using Field = sd::PruningTiming sd::PruningProfile::*;
constexpr std::pair<const char*, Field> fields[] = {
    {"total", &sd::PruningProfile::total},
    {"basis", &sd::PruningProfile::basis},
    {"end0", &sd::PruningProfile::end0},
    {"hom_quotient", &sd::PruningProfile::hom_quotient},
    {"composition", &sd::PruningProfile::composition},
    {"map_lifts", &sd::PruningProfile::map_lifts},
    {"i_phase", &sd::PruningProfile::i_phase},
    {"k_phase", &sd::PruningProfile::k_phase},
    {"i_target_prepare", &sd::PruningProfile::i_target_prepare},
    {"i_shift", &sd::PruningProfile::i_shift},
    {"i_target_reduce", &sd::PruningProfile::i_target_reduce},
    {"i_adaptive_restriction", &sd::PruningProfile::i_adaptive_restriction},
    {"k_adaptive_restriction", &sd::PruningProfile::k_adaptive_restriction},
    {"i_image", &sd::PruningProfile::i_image},
    {"k_image", &sd::PruningProfile::k_image},
    {"i_preimage", &sd::PruningProfile::i_preimage},
    {"i_intersection", &sd::PruningProfile::i_intersection},
    {"k_preimage", &sd::PruningProfile::k_preimage},
    {"i_reduce", &sd::PruningProfile::i_reduce},
    {"k_reduce", &sd::PruningProfile::k_reduce},
    {"i_equal", &sd::PruningProfile::i_equal},
    {"k_equal", &sd::PruningProfile::k_equal},
    {"k_sum", &sd::PruningProfile::k_sum},
    {"k_merge", &sd::PruningProfile::k_merge},
    {"quotient", &sd::PruningProfile::quotient},
    {"presentation", &sd::PruningProfile::presentation},
    {"minimize", &sd::PruningProfile::minimize},
};

double number(const char* text) {
    std::size_t used = 0;
    const std::string value(text);
    const double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result) || result < 0)
        throw std::invalid_argument("Expected a finite nonnegative number");
    return result;
}

double peak_rss_mb() {
    rusage usage_info{};
    if (getrusage(RUSAGE_SELF, &usage_info)) throw std::runtime_error("Cannot read peak RSS");
#ifdef __APPLE__
    return usage_info.ru_maxrss / (1024.0 * 1024.0);
#else
    return usage_info.ru_maxrss / 1024.0;
#endif
}

int main(int argc, char** argv) {
    constexpr const char* usage =
        "Usage: pruning_torus_runtime INPUT.scc EPSILON_PERCENT THREADS MEMORY_MB [--generator-trace] [--standard-i|--parallel-i] [--parallel-k] [--minimize-i-target] [--output-scc=PATH] [--composition-NAME=VALUE ...]\n"
        "Torus runtime trace v1; all-optimizations; hom_algorithm=optimised; threads=0 means auto\n"
        "Generator size trace v1; --standard-i fixes the I target per pass; --parallel-i also processes groups concurrently\n";
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << usage << "Column storage: " << sd::column_storage << '\n';
        return 0;
    }
    std::ostream report(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    report << std::setprecision(17);
    try {
        if (argc < 5) throw std::invalid_argument(usage);
        const double percentage = number(argv[2]), threads = number(argv[3]);
        if (threads != std::trunc(threads) || threads > std::numeric_limits<int>::max())
            throw std::invalid_argument("THREADS must be a nonnegative integer");
        sd::ProcessMemoryGuard memory(number(argv[4]));
        sd::PruningOptions options;
        options.update_i_target = true; // Same all-optimizations configuration as benchmark_runtime.py.
        bool generator_trace = false;
        std::string output_path;
        for (int i = 5; i < argc; ++i) {
            const std::string flag(argv[i]);
            if (flag == "--generator-trace") generator_trace = true;
            else if (flag == "--standard-i") options.update_i_target = false;
            else if (flag == "--parallel-i") options.parallel_i = true;
            else if (flag == "--parallel-k") options.parallel_k = true;
            else if (flag == "--minimize-i-target") options.minimize_i_target = true;
            else if (flag.rfind("--output-scc=", 0) == 0) output_path = flag.substr(13);
            else {
                const auto equal = flag.find('=');
                if (equal == std::string::npos) throw std::invalid_argument("Unknown option: " + flag);
                const auto name = flag.substr(0, equal), value = flag.substr(equal + 1);
                if (name == "--composition-orbit-seed-mode") {
                    if (value != "sum" && value != "random" && value != "hybrid")
                        throw std::invalid_argument("Invalid composition orbit seed mode");
                    options.composition_budget.orbit_seed_mode = value;
                } else if (!sd::parse_composition_count(name, value, options.composition_budget))
                    throw std::invalid_argument("Unknown option: " + flag);
            }
        }
        if (options.parallel_i) options.update_i_target = false;
        if (threads) options.threads = static_cast<int>(threads);
        std::cerr << "Column storage: " << sd::column_storage << "; threads=" << options.threads << '\n';
        sd::print_composition_options(std::cerr, options.composition_budget);
        const auto origin = Clock::now();
        const auto elapsed = [&] { return std::chrono::duration<double>(Clock::now() - origin).count(); };
        const auto load_start = elapsed();
        report << "STEP\tstart\tload_sort\t0\t1\t" << load_start << "\t0" << std::endl;
        sd::Module input(argv[1]);
        input.sort_compatibly();
        const auto epsilon = sd::epsilon_from_percent(input.presentation(), percentage);
        const auto load_end = elapsed();
        report << "STEP\tend\tload_sort\t0\t1\t" << load_end << '\t' << load_end - load_start << std::endl;
        report << "METRICS {\"epsilon_x\":" << epsilon[0] << ",\"epsilon_y\":" << epsilon[1]
               << ",\"threads\":" << options.threads
               << ",\"input_generators\":" << input.num_of_gens()
               << ",\"input_relations\":" << input.num_of_rels()
               << ",\"fast_i\":" << (options.update_i_target ? "true" : "false")
               << ",\"parallel_i\":" << (options.parallel_i ? "true" : "false")
               << ",\"parallel_k\":" << (options.parallel_k ? "true" : "false")
               << ",\"minimize_i_target\":" << (options.minimize_i_target ? "true" : "false") << "}" << std::endl;
        sd::PruningProfile profile;
        std::mutex trace_mutex;
        std::array<std::size_t, std::size(fields)> calls{};
        std::map<std::thread::id, decltype(calls)> active_calls;
        if (generator_trace) profile.trace_generators = [&](const char* object, const char* stage,
                                                           int pass, std::size_t map, std::size_t count) {
            std::lock_guard<std::mutex> lock(trace_mutex);
            report << "GENERATORS\t" << object << '\t' << stage << '\t' << pass << '\t'
                   << map << '\t' << count << '\t' << elapsed() << std::endl;
        };
        profile.trace_step = [&](Field field, bool completed, double seconds) {
            for (std::size_t index = 0; index < std::size(fields); ++index) if (fields[index].second == field) {
                const auto& entry = fields[index];
                const std::string name(entry.first);
                const int pass = name == "i_phase" || name == "k_phase" ? 0 :
                    name.rfind("i_", 0) == 0 ? profile.i_iterations :
                    name.rfind("k_", 0) == 0 ? profile.k_iterations : 0;
                std::lock_guard<std::mutex> lock(trace_mutex);
                auto& call = active_calls[std::this_thread::get_id()][index];
                if (!completed) call = ++calls[index];
                report << "STEP\t" << (completed ? "end" : "start") << '\t' << name << '\t'
                       << pass << '\t' << call
                       << '\t' << elapsed() << '\t' << seconds << std::endl;
                if (completed && field == &sd::PruningProfile::basis)
                    report << "METRICS {\"end_dimension\":" << profile.end_dimension
                           << ",\"quotient_dimension\":" << profile.quotient_dimension
                           << ",\"i_maps\":" << profile.i_maps << ",\"k_maps\":" << profile.k_maps
                           << ",\"peak_rss_mb\":" << peak_rss_mb() << "}" << std::endl;
                return;
            }
        };
        // Match main.cpp's input-copy boundary as well as its actual implementation.
        auto output = sd::pruning(input, epsilon, false, &profile, options);
        report << "RESULT {\"pruning_seconds\":" << profile.total.seconds
               << ",\"epsilon_x\":" << epsilon[0] << ",\"epsilon_y\":" << epsilon[1]
               << ",\"threads\":" << options.threads << ",\"peak_rss_mb\":" << peak_rss_mb()
               << ",\"end_dimension\":" << profile.end_dimension
               << ",\"quotient_dimension\":" << profile.quotient_dimension
               << ",\"i_maps\":" << profile.i_maps << ",\"k_maps\":" << profile.k_maps
               << ",\"i_iterations\":" << profile.i_iterations << ",\"k_iterations\":" << profile.k_iterations
               << ",\"output_generators\":" << output.num_of_gens()
               << ",\"output_relations\":" << output.num_of_rels() << "}" << std::endl;
        if (!output_path.empty()) {
            std::ofstream stream(output_path);
            stream.exceptions(std::ios::failbit | std::ios::badbit);
            output.to_stream(stream);
        }
        sd::print_module_pruning_profile(std::cerr, profile, options);
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
