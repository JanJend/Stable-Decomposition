// Build: cmake --build build --target hom_reduction_benchmark
// Run: build/hom_reduction_benchmark INPUT.scc EPSILON [EPSILON ...] [--percent] [--adaptive]
#include "pruning.hpp"
#include "composition_options.hpp"
#include "process_memory_guard.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace sd = stable_decomposition;
namespace gl = graded_linalg;
using Clock = std::chrono::steady_clock;

template<class Function>
auto timed(double& seconds, Function operation) {
    const auto start = Clock::now();
    auto result = operation();
    seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

void benchmark(std::ostream& out, const std::shared_ptr<const sd::Module>& input,
               const std::vector<sd::Hom>& end,
               double end_seconds, double epsilon, bool adaptive, const gl::CompositionGeneratorOptions& methods) {
    const auto delta = sd::pruning_shift(epsilon);
    const auto eta = sd::Hom::canonical_shift(input, delta);
    out << "\nepsilon=" << epsilon << " shift=" << delta[0] << '\n'
        << "stage\tI_size\tK_size\tseconds\n" << std::flush;
    auto row = [&](const char* stage, std::size_t i, std::size_t k, double seconds) {
        out << stage << '\t' << i << '\t' << k << '\t' << seconds << '\n' << std::flush;
    };

    double hom_seconds, quotient_seconds;
    auto hom = timed(hom_seconds, [&] {
        return gl::module_hom_space_basis(input, eta.target(), gl::HomAlgorithm::hom_exactness);
    });
    row("End_d_basis", hom.size(), hom.size(), hom_seconds);

    std::vector<sd::Hom> transported_end;
    auto q = timed(quotient_seconds, [&] {
        for (const auto& e : end) transported_end.push_back(e.compose(eta));
        return gl::homomorphism_quotient_basis(hom, transported_end);
    });
    row("quotient_basis", q.size(), q.size(), quotient_seconds);

    auto compare_methods = [&](const char* space, const std::vector<sd::Hom>& basis,
                               const std::vector<sd::Hom>& denominator) {
        double setup_seconds;
        auto coordinates = timed(setup_seconds, [&] {
            return basis.empty() || !(methods.sparse_nilpotent || methods.complement_orbit)
                ? std::unique_ptr<gl::PreparedCompositionCoordinates<sd::Mat>>{}
                : std::make_unique<gl::PreparedCompositionCoordinates<sd::Mat>>(basis, denominator);
        });
        row((std::string(space) + "_coordinates").c_str(), basis.size(), basis.size(), setup_seconds);
        for (const auto& method : {std::string("sparse_nilpotent"), std::string("complement_orbit"), std::string("combined")}) {
            const bool sparse = method != "complement_orbit", orbit = method != "sparse_nilpotent";
            if ((sparse && !methods.sparse_nilpotent) || (orbit && !methods.complement_orbit)) continue;
            auto options = methods;
            options.sparse_nilpotent = sparse;
            options.complement_orbit = orbit;
            std::array<gl::CompositionGeneratorStats, 2> stats;
            double seconds;
            const auto selected = timed(seconds, [&] {
                return gl::budgeted_composition_generators(basis, end, denominator, options, &stats, coordinates.get());
            });
            const auto stage = std::string(space) + "_" + method;
            row(stage.c_str(), selected.first.size(), selected.second.size(), seconds);
            for (std::size_t side = 0; side < 2; ++side)
                out << "composition_search stage=" << stage << " side=" << (side ? "K" : "I")
                    << " nilpotent_rank=" << stats[side].nilpotent_rank
                    << " sparse_generators=" << stats[side].sparse_generators
                    << " sparse_pair_visits=" << stats[side].pair_visits
                    << " sparse_products=" << stats[side].products
                    << " orbit_pair_visits=" << stats[side].orbit_pair_visits
                    << " orbit_products=" << stats[side].orbit_products
                    << " orbit_oversized_products=" << stats[side].orbit_oversized_products
                    << " orbit_attempts=" << stats[side].orbit_attempts
                    << " orbit_low_gain_attempts=" << stats[side].orbit_low_gain_attempts
                    << " orbit_seeds=" << stats[side].orbit_seeds
                    << " largest_orbit_seed=" << stats[side].largest_orbit_seed
                    << " rank=" << stats[side].rank << " seeds=" << stats[side].seeds << '\n';
            out << stage << "_route_total_seconds=" << end_seconds + hom_seconds + setup_seconds + seconds
                + (std::string(space) == "quotient" ? quotient_seconds : 0) << '\n' << std::flush;
        }
    };
    compare_methods("quotient", q, transported_end);
    compare_methods("full", hom, {});

    if (adaptive) {
        // Measure the last reduction along real I/K iterates, using the existing
        // profiled algorithm. Its repeated global-basis preparation is separate.
        sd::PruningOptions options;
        options.hom_algorithm = gl::HomAlgorithm::hom_exactness;
        options.composition_generators = true;
        options.composition_budget = methods;
        options.adaptive_restriction = true;
        sd::PruningProfile profile;
        double pair_seconds;
        const auto pair = timed(pair_seconds, [&] {
            return sd::pruning_pair(*input, epsilon, false, &profile, options);
        });
        out << "adaptive_phase\tpasses\tcandidate_visits\tretained_visits\treduction_seconds\n"
            << "I\t" << profile.i_iterations << '\t' << profile.i_restriction_candidates
            << '\t' << profile.i_restriction_kept << '\t' << profile.i_adaptive_restriction.seconds << '\n'
            << "K\t" << profile.k_iterations << '\t' << profile.k_restriction_candidates
            << '\t' << profile.k_restriction_kept << '\t' << profile.k_adaptive_restriction.seconds << '\n'
            << "adaptive_I_target_prepare_seconds=" << profile.i_target_prepare.seconds << '\n'
            << "adaptive_pair_seconds=" << pair_seconds << '\n'
            << "adaptive_recomputed_basis_seconds=" << profile.basis.seconds << '\n'
            << "adaptive_I_generators=" << pair.first.num_of_gens()
            << " adaptive_K_generators=" << pair.second.num_of_gens() << '\n' << std::flush;
    }
}

int main(int argc, char** argv) {
    constexpr const char* usage =
        "Usage: hom_reduction_benchmark INPUT.scc VALUE [VALUE ...] [--percent] [--adaptive]\n"
        "Values are absolute epsilons by default; the shift is 2*epsilon.\n"
        "--percent interprets values as percentages of the largest bounding-box side\n"
        "over all generator and relation degrees: epsilon = VALUE/100 * extent.\n"
        "Values must be finite and nonnegative (1 means 1% with --percent).\n"
        "--memory-limit-mb VALUE optionally limits sampled peak RSS (0 disables).\n"
        "Input is any two-parameter SCC module presentation (e.g. Gaussian-torus H1).\n"
        "Sizes count basis vectors or selected module generators over F2.\n"
        "I uses left End(X)_0 generators; K uses right generators.\n"
        "Compare sparse_nilpotent, complement_orbit, and combined generators.\n"
        "--[no-]sparse-nilpotent-generators and --[no-]complement-orbit-generators select methods.\n"
        "Generator counts are certified but need not be minimum.\n"
        "Method times exclude the separately reported shared coordinate setup.\n"
        "Stage times are incremental; totals include the shared End_0 computation.\n"
        "--adaptive also runs pruning and reports map visits summed over all passes.\n";
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << usage << "--composition-orbit-seed-mode sum|random|hybrid (default: hybrid)\n";
        for (const auto& option : sd::composition_count_options)
            std::cout << "--" << option.name << " N: " << option.help << '\n';
        return 0;
    }
    // Keep library diagnostics/progress on stderr and benchmark tables on stdout.
    std::ostream report(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    try {
        if (argc < 3) throw std::invalid_argument(usage);
        bool adaptive = false, percent = false;
        gl::CompositionGeneratorOptions methods;
        double memory_mb = 0;
        std::vector<double> epsilons;
        for (int i = 2; i < argc; ++i) {
            const std::string raw_argument(argv[i]);
            const auto equals = raw_argument.find('=');
            const std::string argument = raw_argument.substr(0, equals);
            if (argument.rfind("--composition-", 0) == 0) {
                const auto value = equals != std::string::npos ? raw_argument.substr(equals + 1) :
                    (++i < argc ? std::string(argv[i]) : throw std::invalid_argument("Missing value for " + argument));
                if (argument == "--composition-orbit-seed-mode") {
                    if (value != "sum" && value != "random" && value != "hybrid")
                        throw std::invalid_argument("Orbit seed mode must be sum, random, or hybrid");
                    methods.orbit_seed_mode = value;
                } else if (!sd::parse_composition_count(argument, value, methods))
                    throw std::invalid_argument("Unknown option: " + argument);
            } else if (argument == "--adaptive") {
                adaptive = true;
            } else if (argument == "--percent") {
                percent = true;
            } else if (argument == "--sparse-nilpotent-generators" || argument == "--no-sparse-nilpotent-generators") {
                methods.sparse_nilpotent = argument == "--sparse-nilpotent-generators";
            } else if (argument == "--complement-orbit-generators" || argument == "--no-complement-orbit-generators") {
                methods.complement_orbit = argument == "--complement-orbit-generators";
            } else if (argument == "--memory-limit-mb") {
                if (++i == argc) throw std::invalid_argument("Missing memory limit");
                const std::string value(argv[i]);
                std::size_t used = 0;
                memory_mb = std::stod(value, &used);
                if (used != value.size() || !std::isfinite(memory_mb) || memory_mb < 0)
                    throw std::invalid_argument("Memory limit must be finite and nonnegative");
            } else {
                std::size_t used = 0;
                const double epsilon = std::stod(argument, &used);
                if (used != argument.size()) throw std::invalid_argument("Invalid epsilon: " + argument);
                (void)sd::pruning_shift(epsilon);
                epsilons.push_back(epsilon);
            }
        }
        if (epsilons.empty()) throw std::invalid_argument(usage);
        sd::ProcessMemoryGuard memory(memory_mb);
        auto module = std::make_shared<sd::Module>(std::string(argv[1]));
        module->sort_compatibly();
        auto source = module->presentation();
        source.compute_rows_forward();
        double extent = 0;
        if (!source.row_degrees.empty() || !source.col_degrees.empty()) {
            const auto bounds = source.bounding_box();
            extent = std::max(bounds.second[0] - bounds.first[0], bounds.second[1] - bounds.first[1]);
        }
        if (!std::isfinite(extent)) throw std::invalid_argument("Bounding-box extent must be finite");
        report << std::setprecision(17)
               << "input=" << argv[1] << " generators=" << module->num_of_gens()
               << " relations=" << module->num_of_rels() << " storage=" << sd::column_storage << '\n'
               << "bounding_box_extent=" << extent << " epsilon_units=" << (percent ? "percent" : "absolute")
               << '\n' << std::flush;
        sd::print_composition_options(report, methods);
        double end_seconds;
        const auto end = timed(end_seconds, [&] {
            return gl::module_endomorphism_basis(
                std::shared_ptr<const sd::Module>(module), gl::HomAlgorithm::hom_exactness);
        });
        report << "End_0_dimension=" << end.size() << " seconds=" << end_seconds << '\n'
               << "# Incremental seconds; I/K sizes are left/right counts. Generators are heuristic.\n"
               << "# Input loading, sorting and benchmark-only copies are excluded.\n" << std::flush;
        for (const auto value : epsilons) {
            if (percent) report << "\nepsilon_percent=" << value << '\n' << std::flush;
            benchmark(report, module, end, end_seconds, percent ? value * (extent / 100) : value, adaptive, methods);
        }
    } catch (const std::exception& error) {
        std::cerr << "Benchmark error: " << error.what() << '\n';
        return 1;
    }
}
