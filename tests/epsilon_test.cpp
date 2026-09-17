#include <doctest/doctest.h>
#include "epsilon.hpp"
#include "pruning.hpp"
#include "utils.hpp"
#include <limits>

TEST_CASE("epsilon uses all presentation degrees and explicit override") {
    using namespace stable_decomposition;
    Mat empty(0, 0, {}, {}, {});
    CHECK_FALSE(calculate_epsilon_from_matrix(empty).has_value());
    CHECK(get_epsilon(std::nullopt, empty) == doctest::Approx(0.01));
    Mat point(0, 2, {}, {}, {{3,4}, {3,4}});
    CHECK_FALSE(calculate_epsilon_from_matrix(point).has_value());
    Mat free(0, 2, {}, {}, {{-2,1}, {8,4}});
    CHECK(calculate_epsilon_from_matrix(free).value() == doctest::Approx(0.1));
    Mat presentation(1, 1, {{0}}, {{2,20}}, {{-3,0}});
    CHECK(calculate_epsilon_from_matrix(presentation).value() == doctest::Approx(0.2));
    CHECK(get_epsilon(std::nullopt, presentation) == doctest::Approx(0.2));
    CHECK(get_epsilon(0.7, presentation) == doctest::Approx(0.7));
    CHECK(get_epsilon(0.0, presentation) == 0.0);
    Mat cols_only(1, 0, {{}}, {{3,4}}, {});
    CHECK_FALSE(calculate_epsilon_from_matrix(cols_only).has_value());
}

TEST_CASE("epsilon controls the doubled iteration shift and final half shift") {
    using namespace stable_decomposition;
    CHECK(pruning_shift(0.25) == graded_linalg::r2degree(0.5, 0.5));
    // Free generators born at 0 and 1. The map early->late becomes available
    // exactly when 2*epsilon reaches 1, not when epsilon reaches 1.
    Mat births(0, 2, {}, {}, {{0,0}, {1,1}});
    CHECK(End_2epsilon_0(births, 0.25).empty());
    auto additional = End_2epsilon_0(births, 0.5);
    REQUIRE(additional.size() == 1);
    CHECK(additional[0].data == graded_linalg::array<int>({{1}, {}}));
    CHECK(End_2d_0(births, 0.5)[0].data == additional[0].data);
    for (bool quick : {false, true}) {
        // Below the threshold: I is whole, K=0, and final translation is +epsilon.
        auto low = pruning(Module(births), 0.25, quick);
        CHECK(low.presentation().row_degrees == graded_linalg::vec<graded_linalg::r2degree>({{0.25,0.25}, {1.25,1.25}}));
        // At the threshold, invariance under the new map delays the first
        // generator to 1. Both surviving generators then translate to 1.5.
        auto pair = pruning_pair(births, 0.5, quick);
        CHECK(pair.second.get_num_cols() == 0);
        Mat expected_I(2, 2, {{0}, {1}}, {{1,1}, {1,1}}, births.row_degrees);
        CHECK(graded_linalg::present_same_submodule(births, pair.first, expected_I));
        auto high = pruning(Module(births), 0.5, quick);
        CHECK(high.presentation().row_degrees == graded_linalg::vec<graded_linalg::r2degree>({{1.5,1.5}, {1.5,1.5}}));
        CHECK(high.dimension_at({1.25,1.25}) == 0);
        CHECK(high.dimension_at({1.5,1.5}) == 2);
    }
}

TEST_CASE("epsilon validation and historical extractor aliases") {
    using namespace stable_decomposition;
    Mat free(0, 1, {}, {}, {{0,0}});
    for (double value : {-0.1, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::max()}) {
        CHECK_THROWS_AS(pruning_shift(value), std::invalid_argument);
        CHECK_THROWS_AS(get_epsilon(value, free), std::invalid_argument);
        CHECK_THROWS_AS(pruning_pair(free, value), std::invalid_argument);
    }
    CHECK(get_delta(0.25, free) == get_epsilon(0.25, free));
    CHECK(calculate_delta_from_matrix(free) == calculate_epsilon_from_matrix(free));
    CHECK(generate_output_path("dir/input.scc", 0.25) == "dir/input_pru0.2500.scc");
}

TEST_CASE("CLI stores epsilon with a historical delta spelling") {
    auto parse = [](std::vector<std::string> args) {
        std::vector<char*> argv;
        for (auto& arg : args) argv.push_back(arg.data());
        return parse_arguments(static_cast<int>(argv.size()), argv.data());
    };
    CHECK_FALSE(parse({"pruning", "input.scc"}).epsilon.has_value());
    CHECK(parse({"pruning", "input.scc", "--epsilon", "0.125"}).epsilon.value() == 0.125);
    CHECK(parse({"pruning", "input.scc", "--delta", "0.125"}).epsilon.value() == 0.125);
    CHECK_THROWS(parse({"pruning", "input.scc", "--epsilon"}));
    CHECK_THROWS(parse({"pruning", "input.scc", "--epsilon", "0.2junk"}));
    auto flexible = parse({"pruning", "--epsilon=0.25", "--hilbert", "input.scc", "--aida", "-o", "out.scc"});
    CHECK(flexible.epsilon.value() == 0.25);
    CHECK(flexible.hilbert);
    CHECK(flexible.aida);
    CHECK(flexible.output_file == "out.scc");
    CHECK_THROWS(parse({"pruning", "input.scc", "--unknown"}));
    CHECK_THROWS(parse({"pruning", "input.scc", "--image-size", "0"}));
    CHECK(parse({"pruning", "input.scc", "--old"}).old);
    CHECK(parse({"pruning", "input.scc", "--compare"}).compare);
    CHECK(parse({"pruning", "input.scc", "--old", "--quick"}).quick);
    CHECK_THROWS(parse({"pruning", "input.scc", "--old", "--compare"}));
    CHECK_THROWS(parse({"pruning", "input.scc", "--quick"}));
}

TEST_CASE("pruning profiles count actual passes and operations without changing results") {
    using namespace stable_decomposition;
    Mat births(0, 2, {}, {}, {{0,0}, {1,1}});
    for (double epsilon : {0.0, 0.5}) {
        PruningProfile old_profile, new_profile;
        auto matrix_input = births;
        Mat old_result = pruning_profiled(matrix_input, epsilon, false, &old_profile);
        Module new_result = pruning(Module(births), epsilon, false, &new_profile);
        auto plain_input = births;
        Mat old_plain = pruning(plain_input, epsilon, false);
        Module new_plain = pruning(Module(births), epsilon, false);
        CHECK(old_result.data == old_plain.data);
        CHECK(old_result.row_degrees == old_plain.row_degrees);
        CHECK(old_result.col_degrees == old_plain.col_degrees);
        CHECK(new_result.presentation().data == new_plain.presentation().data);
        CHECK(new_result.presentation().row_degrees == new_plain.presentation().row_degrees);
        CHECK(new_result.presentation().col_degrees == new_plain.presentation().col_degrees);
        const int i_passes = epsilon == 0 ? 0 : 2;
        const int k_passes = epsilon == 0 ? 0 : 1;
        for (const auto* p : {&old_profile, &new_profile}) {
            CHECK(p->total.calls == 1);
            CHECK(p->basis.calls == 1);
            CHECK(p->presentation.calls == 1);
            CHECK(p->minimize.calls == 1);
            CHECK(p->i_iterations == i_passes);
            CHECK(p->k_iterations == k_passes);
            CHECK(p->i_equal.calls == i_passes);
            CHECK(p->k_equal.calls == k_passes);
            // This fixture has exactly one additional map when epsilon=0.5.
            CHECK(p->i_preimage.calls == i_passes);
            CHECK(p->k_preimage.calls == k_passes);
            CHECK(p->total.seconds >= p->i_phase.seconds + p->k_phase.seconds);
        }
        CHECK(old_profile.i_intersection.calls == 0);
        CHECK(new_profile.i_intersection.calls == i_passes);
    }
}

TEST_CASE("instrumented matrix copy agrees with original on a presented module") {
    using namespace stable_decomposition;
    for (double epsilon : {0.01, 0.5}) {
        Mat original_input("test1.scc");
        Mat profiled_input = original_input;
        PruningProfile profile;
        Mat original = pruning(original_input, epsilon, false);
        Mat measured = pruning_profiled(profiled_input, epsilon, false, &profile);
        CHECK(measured.data == original.data);
        CHECK(measured.row_degrees == original.row_degrees);
        CHECK(measured.col_degrees == original.col_degrees);
    }
}
