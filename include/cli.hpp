#pragma once
#include "types.hpp"
#include <filesystem>

namespace stable_decomposition {
void write_module(const Module& module, const std::filesystem::path& path);
Module compare_pruning(const Module& input, double epsilon, const std::filesystem::path& prefix);
void write_hilbert_images(Module& input, Module& output, const std::filesystem::path& prefix, int size);
#ifdef PRUNING_WITH_AIDA
void compare_decompositions(const Module& input, const Module& output, const std::filesystem::path& prefix);
#endif
}
