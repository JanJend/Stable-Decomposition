#include "utils.hpp"
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
    app.add_flag("--compare", opts.compare, "Compare matrix/module pruning timings and output invariants (including AIDA when built in); save the module result")
        ->excludes(old);
    app.add_flag("--quick", opts.quick, "Use the original matrix algorithm's quick mode (requires --old)")->needs(old);
    app.add_option("--image-size", opts.image_size, "Hilbert heatmap size in pixels (default: 500)")
        ->check(CLI::Range(128, 2048));
    app.add_flag("--no-output", opts.no_output, "Skip the pruned SCC (optional images/comparison are still saved)");
    app.footer("Options may appear before or after the input; --option=value is supported.\n"
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
    return opts;
}

std::string generate_output_path(const std::string& input, double epsilon) {
    std::filesystem::path p(input);
    std::ostringstream suffix;
    suffix << "_pru" << std::fixed << std::setprecision(4) << epsilon;
    
    std::string new_name = p.stem().string() + suffix.str() + ".scc";
    return (p.parent_path() / new_name).string();
}
