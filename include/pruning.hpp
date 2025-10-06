/**
 * @file pruning.hpp
 * @author Havard Bjerkevik, Jan Jendrysiak, and Fabian Lenzen
 * @brief
 * @version 0.2
 * @date 2025-10-6
 *
 * @copyright
 *
 */

 #pragma once
#ifndef PRUNING_HPP
#define PRUNING_HPP


#include <grlina/graded_linalg.hpp>
#include <vector>

namespace stable_decomposition {
using index_t = int; // Change to large enough int type.

using Mat = graded_linalg::R2GradedSparseMatrix<index_t>;
using namespace graded_linalg;

// Function declarations
void matrix_reduction(vec<Mat>& A, vec<Mat>& B);

vec<Mat> homSpace(Mat& A, Mat& B);

vec<Mat> End_2d_0(Mat &M, double delta);

Mat zero_submodule(const Mat &m);

Mat all_submodule(const Mat &m);

Mat submodule_sum(Mat l, Mat r);

Mat reduce_submodule(Mat M, Mat S);

Mat image(const Mat &f, const Mat &A, const Mat &B, const Mat &U);

Mat shifting_morphism(Mat A, double delta);

bool image_contained_in_image(const Mat &M, const Mat &N);

bool present_same_submodule(const Mat &M, const Mat &A, const Mat &B);

void print_progress(int iteration_I, size_t current, size_t total);

std::vector<Mat> pruning_pair(Mat &M, const double delta);

Mat pruning(Mat &M, const double delta);

// Delta calculation
std::optional<double> calculate_delta_from_matrix(const Mat& M);
double get_delta(std::optional<double> user_delta, const Mat& M);

} // namespace stable_decomposition
#endif // PRUNING_HPP