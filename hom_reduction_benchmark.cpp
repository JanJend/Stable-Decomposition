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
using Coefficients = gl::SparseMatrix<std::ptrdiff_t, gl::DefaultSparseStorage<std::ptrdiff_t>>;
using Clock = std::chrono::steady_clock;

template<class Function>
auto timed(double& seconds, Function operation) {
    const auto start = Clock::now();
    auto result = operation();
    seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

void benchmark(std::ostream& out, const std::shared_ptr<const sd::Module>& input,
               const Coefficients& end,
               double end_seconds, const sd::Epsilon& epsilon, bool adaptive, bool compare, bool dimensions_only,
               bool generators_only, bool time_steps,
               const gl::CompositionGeneratorOptions& methods) {
    const auto delta = sd::pruning_shift(epsilon);
    const auto eta = sd::Hom::canonical_shift(input, delta);
    out << "\nepsilon=" << sd::format_epsilon(epsilon, 17) << " shift=" << sd::format_epsilon(delta, 17) << '\n'
        << "stage\tI_size\tK_size\tseconds\n" << std::flush;
    auto row = [&](const char* stage, std::size_t i, std::size_t k, double seconds) {
        out << stage << '\t' << i << '\t' << k << '\t' << seconds << '\n' << std::flush;
    };

    auto source = input->presentation();
    source.compute_rows_forward();
    double hom_seconds, quotient_seconds;
    auto hom = timed(hom_seconds, [&] {
        return gl::hom_space_coefficients(source, eta.target()->presentation(), sd::PruningOptions{}.hom_algorithm);
    });
    row("End_d_basis", hom.get_num_cols(), hom.get_num_cols(), hom_seconds);

    Coefficients transported_end(0, end.get_num_rows());
    auto q = timed(quotient_seconds, [&] {
        transported_end = end; // Sparse coefficient columns, not graded matrices.
        gl::MapCoefficientContext<sd::Mat>(source, eta.target()->presentation()).normalize(transported_end);
        auto quotient = compare ? hom : std::move(hom);
        gl::reduce_quotient_in_place(transported_end, quotient);
        return quotient;
    });
    row("quotient_basis", q.get_num_cols(), q.get_num_cols(), quotient_seconds);
    if (dimensions_only) return;

    if (!compare) {
        // The production funnel: one quotient, both composition methods in
        // sequence, then restrictions along the actual I and K iterates.
        // Reuse these bases and families; do not prepare Hom a second time.
        double selection_seconds;
        // Match production's lazy coordinate preparation. Keep the legacy
        // setup row; any preparation actually needed is timed in selection.
        row("quotient_coordinates", q.get_num_cols(), q.get_num_cols(), 0);
        if (generators_only) {
            const auto selected = timed(selection_seconds, [&] {
                return gl::budgeted_composition_generators(q, end, transported_end,
                    source, eta.target()->presentation(), methods);
            });
            row("composition_generators", selected.first.get_num_cols(), selected.second.get_num_cols(), selection_seconds);
            return;
        }
        auto families = timed(selection_seconds, [&] {
            auto selected = gl::budgeted_composition_generators(q, end, transported_end,
                source, eta.target()->presentation(), methods);
            return std::make_pair(sd::PruningMapFamily<sd::Mat>(std::move(selected.first), input, eta.target()),
                                  sd::PruningMapFamily<sd::Mat>(std::move(selected.second), input, eta.target()));
        });
        row("composition_generators", families.first.size(), families.second.size(), selection_seconds);
        // Later steps need only eta and the selected families. Release the
        // full coordinate systems before running repeated submodule operations.
        hom = Coefficients(0, end.get_num_rows());
        q = Coefficients(0, end.get_num_rows());
        transported_end = Coefficients(0, end.get_num_rows());
        sd::PruningOptions options;
        options.composition_budget = methods;
        options.adaptive_restriction = true;
        sd::PruningProfile profile;
        auto summary = [&](const char* phase, int passes, std::size_t maximum,
                           std::size_t total, double phase_seconds, double restriction_seconds) {
            out << "restriction_summary phase=" << phase << " passes=" << passes
                << " max_retained=" << maximum << " total_retained=" << total
                << " phase_seconds=" << phase_seconds << " selection_seconds=" << restriction_seconds
                << '\n' << std::flush;
        };
        if (!families.first.empty() || !families.second.empty()) {
            {
                const auto i = sd::compute_pruning_i(eta, families.first, &profile, options, false, &out, time_steps);
            }
            summary("I", profile.i_iterations, profile.i_restriction_max, profile.i_restriction_kept,
                    profile.i_phase.seconds, profile.i_adaptive_restriction.seconds);
            {
                const auto k = sd::compute_pruning_k(eta, families.second, &profile, options, false, &out, time_steps);
            }
        } else {
            // Same empty-family shortcut as pruning_pair: I=X, K=ker(eta).
            out << "empty_quotient_shortcut=1\n" << std::flush;
            summary("I", 0, 0, 0, 0, 0);
        }
        summary("K", profile.k_iterations, profile.k_restriction_max, profile.k_restriction_kept,
                profile.k_phase.seconds, profile.k_adaptive_restriction.seconds);
        return;
    }

    auto compare_methods = [&](const char* space, const Coefficients& basis,
                               const Coefficients& denominator) {
        double setup_seconds;
        auto coordinates = timed(setup_seconds, [&] {
            return basis.get_num_cols() == 0 || !(methods.nilpotent_reduction || methods.complement_orbit)
                ? std::unique_ptr<gl::PreparedCompositionCoordinates<sd::Mat>>{}
                : std::make_unique<gl::PreparedCompositionCoordinates<sd::Mat>>(basis, denominator, source, eta.target()->presentation());
        });
        row((std::string(space) + "_coordinates").c_str(), basis.get_num_cols(), basis.get_num_cols(), setup_seconds);
        for (const auto& method : {std::string("nilpotent_reduction"), std::string("complement_orbit"), std::string("combined")}) {
            const bool nilpotent = method != "complement_orbit", orbit = method != "nilpotent_reduction";
            if ((nilpotent && !methods.nilpotent_reduction) || (orbit && !methods.complement_orbit)) continue;
            auto options = methods;
            options.nilpotent_reduction = nilpotent;
            options.complement_orbit = orbit;
            std::array<gl::CompositionGeneratorStats, 2> stats;
            double seconds;
            const auto selected = timed(seconds, [&] {
                return gl::budgeted_composition_generators(basis, end, denominator, source, eta.target()->presentation(), options, &stats, coordinates.get());
            });
            const auto stage = std::string(space) + "_" + method;
            row(stage.c_str(), selected.first.get_num_cols(), selected.second.get_num_cols(), seconds);
            for (std::size_t side = 0; side < 2; ++side)
                out << "composition_search stage=" << stage << " side=" << (side ? "K" : "I")
                    << " nilpotent_rank=" << stats[side].nilpotent_rank
                    << " certificate_visits=" << stats[side].certificate_visits
                    << " elementary_certificates=" << stats[side].elementary_certificates
                    << " coordinate_preparations=" << stats[side].coordinate_preparations
                    << " nilpotent_generators=" << stats[side].nilpotent_generators
                    << " nilpotent_pair_visits=" << stats[side].pair_visits
                    << " nilpotent_products=" << stats[side].products
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
    compare_methods("full", hom, Coefficients(0, hom.get_num_rows()));

    if (adaptive) {
        // Measure the last reduction along real I/K iterates, using the existing
        // profiled algorithm. Its repeated global-basis preparation is separate.
        sd::PruningOptions options;
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
        "Usage: hom_reduction_benchmark INPUT.scc VALUE [VALUE ...] [--percent] [--dimensions-only|--generators-only] [--compare-methods] [--adaptive]\n"
        "Values are absolute epsilons by default; the shift is 2*epsilon.\n"
        "--percent scales each bounding-box side over all generator and relation degrees.\n"
        "Values must be finite and nonnegative (1 means 1% with --percent).\n"
        "--memory-limit-mb VALUE optionally limits sampled peak RSS (0 disables).\n"
        "Input is any two-parameter SCC module presentation (e.g. Gaussian-torus H1).\n"
        "--percent scales each bounding-box side separately.\n"
        "Sizes count basis vectors or selected module generators over F2.\n"
        "End(X)_0 acts on Q=End(X)_(2 epsilon)/(eta End(X)_0).\n"
        "Default: both composition methods on Q, followed by actual I/K restriction passes.\n"
        "--dimensions-only stops after dim End_0, dim End_(2 epsilon), and dim Q.\n"
        "--generators-only additionally computes left/right E-generators of Q, then stops before pruning.\n"
        "--time-steps appends reduction seconds to each I/K restriction pass (default: off).\n"
        "B_l and B_r count left/right generators of Q, not dimensions of another quotient.\n"
        "--compare-methods compares individual methods on full Hom and Q instead.\n"
        "--[no-]nilpotent-reduction and --[no-]complement-orbit-generators select methods.\n"
        "Generator counts are certified but need not be minimum.\n"
        "Method times exclude the separately reported shared coordinate setup.\n"
        "Stage times are incremental; totals include the shared End_0 computation.\n"
        "--adaptive additionally runs pruning with --compare-methods.\n";
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << usage << "hom_algorithm=" << gl::hom_algorithm_name(sd::PruningOptions{}.hom_algorithm)
                  << " (pruning default)\n"
                  << "end_size_pipeline=1\n"
                  << "--composition-orbit-seed-mode sum|random|hybrid (default: hybrid)\n";
        for (const auto& option : sd::composition_count_options)
            std::cout << "--" << option.name << " N: " << option.help << '\n';
        return 0;
    }
    // Keep library diagnostics/progress on stderr and benchmark tables on stdout.
    std::ostream report(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    try {
        if (argc < 3) throw std::invalid_argument(usage);
        bool adaptive = false, percent = false, compare = false, dimensions_only = false, generators_only = false;
        bool time_steps = false;
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
            } else if (argument == "--compare-methods") {
                compare = true;
            } else if (argument == "--dimensions-only") {
                dimensions_only = true;
            } else if (argument == "--generators-only") {
                generators_only = true;
            } else if (argument == "--time-steps") {
                time_steps = true;
            } else if (argument == "--adaptive") {
                adaptive = true;
            } else if (argument == "--percent") {
                percent = true;
            } else if (argument == "--nilpotent-reduction" || argument == "--no-nilpotent-reduction") {
                methods.nilpotent_reduction = argument == "--nilpotent-reduction";
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
        if (generators_only && (dimensions_only || compare || adaptive))
            throw std::invalid_argument("--generators-only cannot be combined with --dimensions-only, --compare-methods or --adaptive");
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
            return gl::hom_space_coefficients(source, source, sd::PruningOptions{}.hom_algorithm);
        });
        report << "hom_algorithm=" << gl::hom_algorithm_name(sd::PruningOptions{}.hom_algorithm) << '\n'
               << "End_0_dimension=" << end.get_num_cols() << " seconds=" << end_seconds << '\n'
               << "# Incremental seconds; I/K sizes are left/right counts. Generators are heuristic.\n"
               << "# Input loading, sorting and benchmark-only copies are excluded.\n" << std::flush;
        for (const auto value : epsilons) {
            if (percent) report << "\nepsilon_percent=" << value << '\n' << std::flush;
            const auto epsilon = percent ? sd::epsilon_from_percent(source, value) : sd::Epsilon{value, value};
            benchmark(report, module, end, end_seconds, epsilon,
                      adaptive, compare, dimensions_only, generators_only, time_steps, methods);
        }
    } catch (const std::exception& error) {
        std::cerr << "Benchmark error: " << error.what() << '\n';
        return 1;
    }
}
