/** Compatibility names; implementations delegate to Persistence-Algebra. */
#pragma once
#include "types.hpp"

namespace stable_decomposition {
void matrix_reduction(vec<Mat>& A, vec<Mat>& B);
vec<Mat> homSpace(Mat& A, Mat& B);
// Historical name: the actual coordinate shift is (delta, delta).
vec<Mat> End_2d_0(Mat& M, double delta);
Mat zero_submodule(const Mat& M);
Mat all_submodule(const Mat& M);
Mat submodule_sum(Mat A, Mat B);
Mat reduce_submodule(Mat& M, Mat& S);
Mat shifting_morphism(Mat A, double delta);
bool image_contained_in_image(const Mat& A, const Mat& B);
bool present_same_submodule(const Mat& M, const Mat& A, const Mat& B);
/** Image of U in A under the trusted generator lift f: A -> B. */
Mat image(const Mat& f, const Mat& A, const Mat& B, const Mat& U);
} // namespace stable_decomposition
