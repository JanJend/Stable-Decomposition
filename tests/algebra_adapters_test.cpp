#include <doctest/doctest.h>
#include "pruning.hpp"

TEST_CASE("legacy algebra adapters keep their public names") {
    namespace sd = stable_decomposition;
    using sd::Mat;
    Mat free(0, 1, {}, {}, {{0,0}});
    Mat whole(1, 1, {{0}}, {{0,0}}, {{0,0}});
    auto zero = sd::zero_submodule(free);
    CHECK(zero.get_num_cols() == 0);
    CHECK(sd::all_submodule(free).data == whole.data);
    CHECK(sd::submodule_sum(whole, whole).get_num_cols() == 1);
    CHECK(sd::reduce_submodule(free, whole).data == whole.data);
    CHECK(sd::image_contained_in_image(zero, whole));
    CHECK_FALSE(sd::image_contained_in_image(whole, zero));
    CHECK(sd::present_same_submodule(free, whole, sd::all_submodule(free)));
    auto shift = sd::shifting_morphism(free, 1);
    CHECK(shift.row_degrees == graded_linalg::vec<graded_linalg::r2degree>({{-1,-1}}));
    CHECK(shift.col_degrees == free.row_degrees);
    CHECK(sd::homSpace(free, free).size() == 1);
    CHECK(free.rows_computed);
    CHECK(sd::image(whole, free, free, whole).data == whole.data);
    sd::vec<Mat> A{whole}, B{whole};
    sd::matrix_reduction(A, B);
    CHECK(B.empty());
}

TEST_CASE("empty additional lift family terminates quick pruning") {
    namespace sd = stable_decomposition;
    sd::Mat square(2, 1, {{0}, {0}}, {{0,1}, {1,0}}, {{0,0}});
    CHECK(sd::End_2d_0(square, 0).empty());
    auto pair = sd::pruning_pair(square, 0, true);
    CHECK(pair.first.data == graded_linalg::array<int>({{0}}));
    CHECK(pair.second.get_num_cols() == 0);
    auto result = sd::pruning(sd::PModule(square), 0, true);
    CHECK(result.dimension_at({0.5,0.5}) == 1);
    CHECK(result.dimension_at({1,0}) == 0);
}

TEST_CASE("timing utility joins on success and exception") {
    CHECK(stable_decomposition::timed_with_progress("value", [] { return 7; }) == 7);
    int value = 0;
    stable_decomposition::timed_with_progress("void", [&] { value = 3; });
    CHECK(value == 3);
    auto owned = stable_decomposition::timed_with_progress("move-only", [] {
        return std::make_unique<int>(8);
    });
    CHECK(*owned == 8);
    int& reference = stable_decomposition::timed_with_progress("reference", [&]() -> int& { return value; });
    CHECK(&reference == &value);
    CHECK_THROWS_AS(stable_decomposition::timed_with_progress("failure", []() -> int {
        throw std::runtime_error("expected");
    }), std::runtime_error);
    stable_decomposition::print_progress(1, 0, 0);
    stable_decomposition::print_progress(1, 1, 1);
}
