#include "utils.hpp"
#include "types.hpp"
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
    app.add_flag("--fuse-i-preimage", opts.pruning_options.fuse_i_preimage,
                 "Module: combine the I preimage and intersection")->excludes(old)->excludes(compare);
    app.add_flag("--defer-k-image-minimization", opts.pruning_options.defer_k_image_minimization,
                 "Module: minimize K images only when summing")->excludes(old)->excludes(compare);
    app.add_flag("--projected-preimages", opts.pruning_options.projected_preimages,
                 "Module: track only needed kernel coordinates in two-parameter preimages")->excludes(old)->excludes(compare);
    app.add_flag("--dense-kernel-pivots", opts.pruning_options.dense_kernel_pivots,
                 "Module: use indexed pivot lookup in two-parameter iteration kernels")->excludes(old)->excludes(compare);
    app.add_flag("--throttle-progress", opts.pruning_options.throttle_progress,
                 "Module: update progress only when its percentage changes")->excludes(old)->excludes(compare);
    app.add_flag("--composition-generators", opts.pruning_options.composition_generators,
                 "Module: use certified post/precomposition generators of the Hom quotient")->excludes(old)->excludes(compare);
    app.add_flag("--amortize-reductions", opts.pruning_options.amortize_reductions,
                 "Module: batch K image sums before reducing generators")->excludes(old)->excludes(compare);
    app.add_flag("--hom-classes", opts.pruning_options.hom_classes,
                 "Module: construct actual Hom quotient classes without composition selection")->excludes(old)->excludes(compare);
    app.add_flag("--phase-quotient", opts.pruning_options.phase_quotient,
                 "Module: remove dependent map restrictions separately in each I/K pass")->excludes(old)->excludes(compare);
    app.add_flag("--reuse-target-reductions", opts.pruning_options.reuse_target_reductions,
                 "Module: reduce the fixed I-pass target generators once for all preimages")->excludes(old)->excludes(compare);
    bool all_optimizations = false;
    app.add_flag("--all-optimizations", all_optimizations, "Module: enable all ten optimization switches")
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
    if (all_optimizations) opts.pruning_options = {true, true, true, true, true, true, true, true, true, true};
    return opts;
}

std::string generate_output_path(const std::string& input, double epsilon) {
    std::filesystem::path p(input);
    std::ostringstream suffix;
    suffix << "_pru" << std::fixed << std::setprecision(4) << epsilon;
    
    std::string new_name = p.stem().string() + suffix.str() + ".scc";
    return (p.parent_path() / new_name).string();
}
