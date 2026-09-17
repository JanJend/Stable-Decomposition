
#include "include/pruning.hpp"
#include "include/utils.hpp"
#include <iostream>
#include <fstream>

using namespace graded_linalg;
using namespace stable_decomposition;

int main(int argc, char** argv) try {
    auto opts = parse_arguments(argc, argv);
    if (opts.input_file.empty()) return 1;
    
    stable_decomposition::Module M(opts.input_file);
    M.sort_compatibly();
    const double epsilon = get_epsilon(opts.epsilon, M.presentation());
    
    std::cout << "Computing pruning of " << opts.input_file 
              << " (epsilon=" << epsilon << ", shift=" << 2 * epsilon << ")" << std::endl;
    
    stable_decomposition::Module Pru_M = pruning(std::move(M), epsilon, true);
    
    if (!opts.no_output) {
        std::string output_path = generate_output_path(opts.input_file, epsilon);
        std::ofstream output_file(output_path);
        
        if (!output_file.is_open()) {
            std::cerr << "Error: Unable to open " << output_path << std::endl;
            return 1;
        }
        
        Pru_M.to_stream(output_file);
        std::cout << "Saved to: " << output_path << std::endl;
    }
    
    return 0;
}
catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << std::endl;
    return 1;
}
