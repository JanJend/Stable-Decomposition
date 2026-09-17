#include "cli.hpp"
#include "pruning.hpp"
#include "utils.hpp"
#include <iostream>

using namespace stable_decomposition;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
    try {
        const auto options = parse_arguments(argc, argv);
        if (options.help) return 0;
#ifndef PRUNING_WITH_AIDA
        if (options.aida)
            throw std::runtime_error("AIDA support is disabled; rebuild with -DPRUNING_WITH_AIDA=ON");
#endif

        // Read the input module and choose the pruning scale.
        Module input(options.input_file);
        input.sort_compatibly();
        const double epsilon = get_epsilon(options.epsilon, input.presentation());

        // Optional images and decompositions share the output file's directory and stem.
        const fs::path output_path = options.output_file.empty()
            ? generate_output_path(options.input_file, epsilon) : options.output_file;
        if (!options.no_output && fs::weakly_canonical(output_path) == fs::weakly_canonical(options.input_file))
            throw std::invalid_argument("Choose an output path different from the input");
        auto output_prefix = output_path;
        output_prefix.replace_extension();
        if ((!options.no_output || options.hilbert || options.aida) && !output_path.parent_path().empty())
            fs::create_directories(output_path.parent_path());

        std::cout << "Computing pruning of " << options.input_file << " (epsilon=" << epsilon << ")\n";
        Module output;
        if (options.compare) {
            output = compare_pruning(input, epsilon, output_prefix);
        } else if (options.old) {
            Mat presentation = input.presentation(); // Matrix pruning modifies its input.
            output = Module(pruning(presentation, epsilon, options.quick));
        } else {
            output = pruning(input, epsilon, false);
        }

        if (!options.no_output) {
            write_module(output, output_path);
            std::cout << "Saved to: " << (options.compare ? output_path.filename() : output_path).string() << '\n';
        }
        if (options.hilbert)
            write_hilbert_images(input, output, output_prefix, options.image_size);
#ifdef PRUNING_WITH_AIDA
        if (options.aida)
            compare_decompositions(input, output, output_prefix);
#endif
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\nRun with --help for usage.\n";
        return 1;
    }
}
