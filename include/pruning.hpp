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

namespace stable_decomposition {
using index_t = int; // Change to large enough int type.

using Mat = graded_linalg::R2GradedSparseMatrix<index_t>;
using graded_linalg::vec;

void matrix_reduction(vec<Mat>& A, vec<Mat>& B);

vec<Mat> homSpace(Mat& A, Mat& B);

/**
 * Computes a basis of M(2δ).
 *
 * @param M presentation matrix for the module
 * @param delta the shift 
 */
vec<Mat> End_2d_0(Mat &M, double delta);

/**
 * Generators matrix for the zero submodule of a module.
 * 
 * @param M m×n presentation matrix of the ambient module.
 * @returns m×0 matrix with the same row grades as m.
 */
Mat zero_submodule(const Mat &M);

/**
 * Generators matrix for a module as submodule of itself.
 * 
 * @param M m×n presentation matrix of the ambient module.
 * @returns m×m identity matrix with the same row and column grades as row grades of M.
 */
Mat all_submodule(const Mat &m);

/**
 * Juxtapose matrices and reduce.
 */
Mat submodule_sum(Mat l, Mat r);

/**
 * Simplifies the presentation of S as a submodule of M.
 * @param S generators of a submodule S ⊆ M
 * @param M presentation of the ambient module M
 */
Mat reduce_submodule(Mat& M, Mat& S);

Mat image(const Mat &f, const Mat &A, const Mat &B, const Mat &U);

/**
 * Matrix representing the canonical morphism M → M(δ).
 * 
 * @param M m× n presentation matrix of the module M, with vector of row grades r 
 * @returns m× m identity matrix, with row grades r(δ) and column grades r.
 */
Mat shifting_morphism(Mat A, double delta);

bool image_contained_in_image(const Mat &M, const Mat &N);

/**
 * True if two matrices generate the same submodule of a given one.
 * 
 * @param M presentation of a module M
 * @param A generators of a submodule A ⊆ M
 * @param B generators of a submodule B ⊆ M
 * @returns true of the submodules A, B ⊆ M are identical.
 */
bool present_same_submodule(const Mat &M, const Mat &A, const Mat &B);

void print_progress(int iteration_I, size_t current, size_t total);

/**
 * Computes the *pruning pair* (I, K) of the module presented by the matrix M with respect to the shift δ.
 * 
 * The pruning pair of M is computed iteratively by
 * - I = ⋂ᵢ Iᵢ for I₀ = M and Iᵢ₊₁ = ⋂_f f⁻¹(sh(Iᵢ, 2δ)),
 * - K = ⋃ᵢ Kᵢ for K₀ = 0 and Kᵢ₊₁ = ∑_f sh⁻¹(f(Kᵢ), 2δ),
 * where f runs over a basis of End(M)_{2δ}.
 * @param M presentation of a module M
 * @returns δ-pruning pair [I, K] of M
 */
std::pair<Mat, Mat> pruning_pair(Mat &M, const double delta);

/**
 * Computes the δ-pruning of a module.
 *
 * The δ-pruning of a module M with δ-pruning pair K ⊆ I ⊆ M is the module I/K.
 * @param M presentation matrix of the module
 * @returns a presentation matrix for I/K
 */
Mat pruning(Mat &M, const double delta);

// Delta calculation
std::optional<double> calculate_delta_from_matrix(const Mat& M);
double get_delta(std::optional<double> user_delta, const Mat& M);

} // namespace stable_decomposition
#endif // PRUNING_HPP