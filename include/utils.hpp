#pragma once
#include <optional>
#include <string>

struct ProgramOptions {
    std::string input_file, output_file;
    std::optional<double> epsilon;
    bool help = false;
    bool no_output = false;
    bool hilbert = false;
    bool aida = false;
    int image_size = 500;
};

ProgramOptions parse_arguments(int argc, char** argv);
std::string generate_output_path(const std::string& input, double epsilon);
