#pragma once

#include <string>
#include <optional>

struct ProgramOptions {
    std::string input_file;
    std::optional<double> delta;  // None if not specified
    bool no_output = false;
    bool no_timers = false;
};

ProgramOptions parse_arguments(int argc, char** argv);
std::string generate_output_path(const std::string& input, double delta);
