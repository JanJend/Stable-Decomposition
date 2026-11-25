// Small subprocess worker. Python supplies independent deadlines for each phase.
#include "pruning.hpp"
#include <grlina/invariants/hilbert_function.hpp>
#include <grlina/module/epsilon.hpp>
#include <aida_interface.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace sd = stable_decomposition;
using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

sd::Epsilon parse_epsilon(const char* text) {
    std::size_t used = 0;
    const std::string value(text);
    const double x = std::stod(value, &used);
    sd::Epsilon epsilon{x, x};
    if (used != value.size()) {
        if (value[used] != ',') throw std::invalid_argument("Invalid epsilon");
        const auto tail = value.substr(used + 1);
        epsilon[1] = std::stod(tail, &used);
        if (used != tail.size()) throw std::invalid_argument("Invalid epsilon");
    }
    (void)sd::pruning_shift(epsilon);
    return epsilon;
}

int main(int argc, char** argv) {
    // Algorithm progress goes to stderr; stdout contains only machine-readable metrics.
    std::ostream report(std::cout.rdbuf());
    std::cout.rdbuf(std::cerr.rdbuf());
    report << std::setprecision(17);
    try {
        if (argc < 3) throw std::invalid_argument(
            "Usage: pruning_experiment inspect INPUT | prune INPUT EPSILON TEMP_OUTPUT | components INPUT [EPSILON ...]");
        const std::string mode(argv[1]);
        if (mode == "inspect" && argc == 3) {
            sd::Module input(argv[2], true);
            auto lower = graded_linalg::r2degree{0, 0}, upper = lower;
            if (input.num_of_gens() || input.num_of_rels())
                std::tie(lower, upper) = input.bounding_box();
            report << "RESULT {\"generators\":" << input.num_of_gens()
                   << ",\"relations\":" << input.num_of_rels()
                   << ",\"bbox_min_x\":" << lower[0] << ",\"bbox_min_y\":" << lower[1]
                   << ",\"bbox_max_x\":" << upper[0] << ",\"bbox_max_y\":" << upper[1]
                   << ",\"bbox_extent\":" << std::max(upper[0] - lower[0], upper[1] - lower[1])
                   << "}\n" << std::flush;
            const auto start = Clock::now();
            input.compute_projective_resolution();
            const auto thickness = input.hilbert_function().maximum;
            report << "RESULT {\"thickness\":" << thickness
                   << ",\"thickness_seconds\":" << seconds(start) << "}\n";
        } else if (mode == "prune" && argc == 5) {
            sd::Module input(argv[2], true);
            const auto epsilon = parse_epsilon(argv[3]);
            sd::PruningProfile profile;
            const auto start = Clock::now();
            const auto result = sd::pruning(std::move(input), epsilon, false, &profile);
            const double elapsed = seconds(start);
            // Only the temporary handoff to AIDA is written, outside the timer.
            std::ofstream output(argv[4]);
            output.exceptions(std::ios::failbit | std::ios::badbit);
            result.presentation().to_stream(output);
            output.close();
            report << "RESULT {\"pruning_seconds\":" << elapsed
                   << ",\"i_iterations\":" << profile.i_iterations
                   << ",\"k_iterations\":" << profile.k_iterations << "}\n";
        } else if (mode == "components") {
            // AIDA's editable coordinate algorithm expects vector columns. This
            // explicit boundary avoids changing the pruning backend's storage.
            using Matrix = graded_linalg::CoordGradedSparseMatrix<double, 2, int, void,
                std::vector<std::vector<int>>>;
            using Module = graded_linalg::Module<Matrix>;
            Module input(argv[2], true);
            input.minimize();
            aida::AIDA_functor<> decomposer;
            aida::Block_list_t<Matrix> blocks;
            const auto start = Clock::now();
            if (input.num_of_gens()) decomposer(input.presentation(), blocks);
            const double elapsed = seconds(start);
            // Classify lifetimes separately from decomposition timing. A component
            // dies at a vector shift iff the corresponding shifted image is zero.
            const auto classification_start = Clock::now();
            std::vector<std::size_t> short_epsilon(argc - 3), short_twice(argc - 3);
            for (const auto& block : blocks) {
                auto summand = std::make_shared<const Module>(static_cast<const Matrix&>(block));
                for (int i = 3; i < argc; ++i) {
                    const auto epsilon = parse_epsilon(argv[i]);
                    auto image = graded_linalg::Submodule<Matrix>::whole(summand);
                    image.shift_generators(epsilon);
                    short_epsilon[i - 3] += image.is_zero();
                    image.shift_generators(epsilon);
                    short_twice[i - 3] += image.is_zero();
                }
            }
            report << "RESULT {\"components\":" << blocks.size()
                   << ",\"aida_seconds\":" << elapsed
                   << ",\"classification_seconds\":" << seconds(classification_start)
                   << ",\"short_epsilon\":[";
            for (std::size_t i = 0; i < short_epsilon.size(); ++i)
                report << (i ? "," : "") << short_epsilon[i];
            report << "],\"short_2epsilon\":[";
            for (std::size_t i = 0; i < short_twice.size(); ++i)
                report << (i ? "," : "") << short_twice[i];
            report << "]}\n";
        } else throw std::invalid_argument("Unknown mode or wrong argument count");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Experiment error: " << error.what() << '\n';
        return 1;
    }
}
