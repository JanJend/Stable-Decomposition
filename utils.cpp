#include "include/utils.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <filesystem>  

ProgramOptions parse_arguments(int argc, char** argv) {
    ProgramOptions opts;
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <file_path> [options]\n"
                  << "Options:\n"
                  << "  --delta <value>   Set delta value (default: 0.02)\n"
                  << "  --no-output       Skip saving output file\n"
                  << "  --no-timers       Disable timing output\n";
        opts.input_file = "/home/wsljan/MP-Workspace/mpm_generation/skyscraper/handcrafted_example/two_circles.scc";
        return opts;
    }
    
    opts.input_file = argv[1];
    
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-output" || arg == "-no-output") {
            opts.no_output = true;
        } else if (arg == "--no-timers" || arg == "-no-timers") {
            opts.no_timers = true;
        } else if ((arg == "--delta" || arg == "-delta") && i + 1 < argc) {
            opts.delta = std::stod(argv[++i]);
        }
    }
    
    return opts;
}

std::string generate_output_path(const std::string& input, double delta) {
    std::filesystem::path p(input);
    std::ostringstream suffix;
    suffix << "_pru" << std::fixed << std::setprecision(4) << delta;
    
    std::string new_name = p.stem().string() + suffix.str() + p.extension().string();
    return (p.parent_path() / new_name).string();
}