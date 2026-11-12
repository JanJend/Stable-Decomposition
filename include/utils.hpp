#pragma once
#include <optional>
#include <string>
#include "pruning_options.hpp"

struct ProgramOptions {
    double memory_limit_mb = 0;
    std::string input_file, output_file;
    std::optional<double> epsilon;
    bool help = false;
    bool no_output = false;
    bool hilbert = false;
    bool aida = false;
    bool old = false;
    bool compare = false;
    bool quick = false;
    bool profile = false;
    stable_decomposition::PruningOptions pruning_options;
    int image_size = 500;
};

ProgramOptions parse_arguments(int argc, char** argv);
std::string generate_output_path(const std::string& input, double epsilon);
