#include "utils.hpp"
#include "types.hpp"
#include "composition_options.hpp"
#include <CLI11.hpp>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <filesystem>  
#include <stdexcept>

ProgramOptions parse_arguments(int argc, char** argv) {
    ProgramOptions opts;
    CLI::App app{"Prune a two-parameter persistence module: SCC input -> SCC output."};
    app.set_version_flag("-v,--version", "Stable-Decomposition " PRUNING_VERSION);
    app.add_option("input", opts.input_file, "Input SCC presentation or projective resolution")->required();
    app.add_option("-e,--epsilon,--delta", opts.epsilon,
                   "Nonnegative epsilon (default: 1% of degree extent; fallback 0.01)");
    app.add_option("-o,--output", opts.output_file, "Output SCC path (default: <input>_pru<epsilon>.scc)");
    app.add_flag("--hilbert", opts.hilbert, "Save input/output Hilbert PNGs with shared axes and colour scale");
    app.add_flag("--aida", opts.aida, "Decompose input/output with AIDA; save decompositions and comparison");
    auto* old = app.add_flag("--old", opts.old, "Use the original matrix pruning implementation");
    auto* compare = app.add_flag("--compare", opts.compare, "Compare matrix/module pruning timings and output invariants (including AIDA when built in); save the module result")
        ->excludes(old);
    app.add_flag("--quick", opts.quick, "Use the original matrix algorithm's quick mode (requires --old)")->needs(old);
    app.add_flag("--profile", opts.profile, "Print module pruning operation timings")
        ->excludes(old)->excludes(compare);
    app.add_flag("--projected-preimages,!--no-projected-preimages", opts.pruning_options.projected_preimages,
                 "Module: track only needed kernel coordinates (default: on)")->excludes(old)->excludes(compare);
    app.add_flag("--composition-generators,!--no-composition-generators", opts.pruning_options.composition_generators,
                 "Module: use certified composition generators of the Hom quotient (default: on)")->excludes(old)->excludes(compare);
    app.add_flag("--sparse-nilpotent-generators,!--no-sparse-nilpotent-generators",
                 opts.pruning_options.composition_budget.sparse_nilpotent,
                 "Module: nilpotent actions and small sparse seeds (default: on)")->excludes(old)->excludes(compare);
    app.add_flag("--complement-orbit-generators,!--no-complement-orbit-generators",
                 opts.pruning_options.composition_budget.complement_orbit,
                 "Module: orbit compression of sums of remaining basis maps (default: on)")->excludes(old)->excludes(compare);
    for (const auto& option : stable_decomposition::composition_count_options)
        app.add_option("--" + std::string(option.name), opts.pruning_options.composition_budget.*(option.member), option.help)
            ->check(CLI::NonNegativeNumber)->excludes(old)->excludes(compare);
    app.add_option("--composition-orbit-seed-mode", opts.pruning_options.composition_budget.orbit_seed_mode,
                   "Orbit candidates: sum, random, hybrid (alternating broad/random; default)")
        ->check(CLI::IsMember({"sum", "random", "hybrid"}))->excludes(old)->excludes(compare);
    app.add_flag("--amortize-reductions,!--no-amortize-reductions", opts.pruning_options.amortize_reductions,
                 "Module: batch K image sums before reducing generators (default: on)")->excludes(old)->excludes(compare);
    app.add_flag("--full-end-space", opts.pruning_options.full_end_space,
                 "Module: use End_(2 epsilon) without quotienting by shifted End_0 (default: off)")->excludes(old)->excludes(compare);
    std::string hom_algorithm = "optimised";
    app.add_option("--hom-algorithm", hom_algorithm,
                   "Module Hom algorithm: optimised (A-1/2, default), hom-exactness (B), full-restriction (A)")
        ->check(CLI::IsMember({"optimised", "full-restriction", "hom-exactness"}))
        ->excludes(old)->excludes(compare);
    app.add_flag("--adaptive-restriction,!--no-adaptive-restriction", opts.pruning_options.adaptive_restriction,
                 "Module: remove dependent restrictions in each I/K pass (default: on)")->excludes(old)->excludes(compare);
    app.add_flag("--reuse-target-reductions,!--no-reuse-target-reductions", opts.pruning_options.reuse_target_reductions,
                 "Module: reduce the fixed I-pass target generators once (default: on)")->excludes(old)->excludes(compare);
    bool all_optimizations = false;
    app.add_flag("--all-optimizations", all_optimizations, "Module: enable all optional optimizations and the Hom quotient (already the default); keep the selected Hom algorithm and search budgets")
        ->excludes(old)->excludes(compare);
    app.add_option("--image-size", opts.image_size, "Hilbert heatmap size in pixels (default: 500)")
        ->check(CLI::Range(128, 2048));
    app.add_flag("--no-output", opts.no_output, "Skip the pruned SCC (optional images/comparison are still saved)");
    app.add_option("--memory-limit-mb", opts.memory_limit_mb,
                   "Sample peak RSS and stop above this many MiB (0 disables the limit)")
        ->check(CLI::Range(0.0, 1e9));
    app.footer(std::string("Column storage: ") + stable_decomposition::column_storage + "\n"
               "Options may appear before or after the input; --option=value is supported.\n"
               "Example: pruning input.scc -e 0.05 --hilbert --aida -o results/pruned.scc");
    try {
        app.parse(argc, argv);
    } catch (const CLI::Success& message) {
        app.exit(message);
        opts.help = true;
        return opts;
    } catch (const CLI::ParseError& error) {
        throw std::invalid_argument(error.what());
    }
    if (opts.epsilon && (!std::isfinite(*opts.epsilon) || *opts.epsilon < 0 || !std::isfinite(2 * *opts.epsilon)))
        throw std::invalid_argument("epsilon and 2*epsilon must be finite and nonnegative");
    using graded_linalg::HomAlgorithm;
    opts.pruning_options.hom_algorithm = hom_algorithm == "full-restriction" ? HomAlgorithm::full_restriction
        : hom_algorithm == "hom-exactness" ? HomAlgorithm::hom_exactness : HomAlgorithm::optimised;
    if (all_optimizations) {
        auto& options = opts.pruning_options;
        options.projected_preimages = options.composition_generators = options.amortize_reductions = true;
        options.adaptive_restriction = options.reuse_target_reductions = true;
        options.composition_budget.sparse_nilpotent = options.composition_budget.complement_orbit = true;
        options.full_end_space = false;
    }
    return opts;
}

std::string generate_output_path(const std::string& input, double epsilon) {
    std::filesystem::path p(input);
    std::ostringstream suffix;
    suffix << "_pru" << std::fixed << std::setprecision(4) << epsilon;
    
    std::string new_name = p.stem().string() + suffix.str() + ".scc";
    return (p.parent_path() / new_name).string();
}
