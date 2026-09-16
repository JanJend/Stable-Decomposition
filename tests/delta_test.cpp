#include <doctest/doctest.h>
#include "delta.hpp"

TEST_CASE("delta uses all presentation degrees and explicit override") {
    using namespace stable_decomposition;
    Mat empty(0, 0, {}, {}, {});
    CHECK_FALSE(calculate_delta_from_matrix(empty).has_value());
    CHECK(get_delta(std::nullopt, empty) == doctest::Approx(0.01));
    Mat point(0, 2, {}, {}, {{3,4}, {3,4}});
    CHECK_FALSE(calculate_delta_from_matrix(point).has_value());
    Mat free(0, 2, {}, {}, {{-2,1}, {8,4}});
    CHECK(calculate_delta_from_matrix(free).value() == doctest::Approx(0.1));
    Mat presentation(1, 1, {{0}}, {{2,20}}, {{-3,0}});
    CHECK(calculate_delta_from_matrix(presentation).value() == doctest::Approx(0.2));
    CHECK(get_delta(std::nullopt, presentation) == doctest::Approx(0.2));
    CHECK(get_delta(0.7, presentation) == doctest::Approx(0.7));
    CHECK(get_delta(0.0, presentation) == 0.0);
    Mat cols_only(1, 0, {{}}, {{3,4}}, {});
    CHECK_FALSE(calculate_delta_from_matrix(cols_only).has_value());
}
