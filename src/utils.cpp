#include "utils.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <filesystem>  
#include <stdexcept>

ProgramOptions parse_arguments(int argc, char** argv) {
    ProgramOptions opts;
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file_path> [options]\n"
                  << "Options:\n"
                  << "  --epsilon <value> Set epsilon (default: 1% of degree extent; fallback 0.01)\n"
                  << "  --delta <value>   Compatibility alias for --epsilon\n"
                  << "  --no-output       Skip saving output file\n"
                  << "  --no-timers       Disable timing output\n";
        return opts;
    }
    
    opts.input_file = argv[1];
    
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-output" || arg == "-no-output") {
            opts.no_output = true;
        } else if (arg == "--no-timers" || arg == "-no-timers") {
            opts.no_timers = true;
        } else if (arg == "--epsilon" || arg == "-epsilon" || arg == "--delta" || arg == "-delta") {
            if (i + 1 >= argc) throw std::invalid_argument("Missing value for epsilon");
            const std::string value = argv[++i];
            std::size_t used = 0;
            opts.epsilon = std::stod(value, &used);
            if (used != value.size()) throw std::invalid_argument("Invalid epsilon value: " + value);
        }
    }
    
    return opts;
}

std::string generate_output_path(const std::string& input, double epsilon) {
    std::filesystem::path p(input);
    std::ostringstream suffix;
    suffix << "_pru" << std::fixed << std::setprecision(4) << epsilon;
    
    std::string new_name = p.stem().string() + suffix.str() + p.extension().string();
    return (p.parent_path() / new_name).string();
}
