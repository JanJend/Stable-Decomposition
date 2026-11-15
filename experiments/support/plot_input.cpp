#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "types.hpp"
#include <grlina/invariants/hilbert_function.hpp>
#include <grlina/invariants/hilbert_plot.hpp>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    std::ostream report(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    try {
        if (argc != 4)
            throw std::invalid_argument("Usage: plot_input INPUT OUTPUT_PNG TITLE");
        const auto start = std::chrono::steady_clock::now();
        stable_decomposition::Module input(argv[1], true);
        const auto degrees = graded_linalg::support_degrees(input);
        double xmax = 0;
        for (const auto& degree : degrees) {
            if (!std::isfinite(degree[0]) || !std::isfinite(degree[1]) || degree[0] < 0)
                throw std::invalid_argument("Expected finite radius/codensity grades");
            xmax = std::max(xmax, degree[0]);
        }
        xmax = xmax > 0 ? xmax * 1.05 : 1;
        const auto xs = graded_linalg::image_detail::uniform_axis(0.0, xmax, 1000);
        const auto ys = graded_linalg::image_detail::uniform_axis(-1.02, 0.02, 700);
        input.compute_projective_resolution();
        // Evaluate exactly at the pixels; avoid a potentially huge Cartesian
        // table on every critical grade of the largest presentations.
        const auto values = graded_linalg::euler_grid_r2(input.proj_res_, xs, ys);
        int maximum = 0;
        for (const auto& column : values) for (const auto value : column) {
            if (value < 0) throw std::runtime_error("Negative Hilbert dimension");
            maximum = std::max(maximum, value);
        }
        const auto png = graded_linalg::image_detail::heatmap_image(
            xs, ys, values, argv[3], "Dimension", true, 1, maximum);
        const std::filesystem::path output(argv[2]);
        std::filesystem::create_directories(output.parent_path());
        graded_linalg::save_png(png, output.string());
        report << std::setprecision(17)
               << "{\"radius_max\":" << xmax
               << ",\"sampled_max_dimension\":" << maximum
               << ",\"far_corner_dimension\":" << values.back().back()
               << ",\"resolution_maps\":" << input.proj_res_.size()
               << ",\"seconds\":" << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start).count() << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
