
#include "include/pruning.hpp"
#include "include/utils.hpp"
#include <iostream>
#include <fstream>

using namespace graded_linalg;
using namespace stable_decomposition;

int main(int argc, char** argv) {
    auto opts = parse_arguments(argc, argv);
    
    Mat M(opts.input_file);
    M.sort_columns_lexicographically();
    M.sort_rows_lexicographically();
    double delta = get_delta(opts.delta, M);
    
    std::cout << "Computing pruning of " << opts.input_file 
              << " (delta=" << delta << ")" << std::endl;
    
    Mat Pru_M = pruning(M, delta);
    
    if (!opts.no_output) {
        std::string output_path = generate_output_path(opts.input_file, delta);
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