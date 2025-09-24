/**
 * @file pruning.hpp
 * @author Jan Jendrysiak and Fabian Lenzen
 * @brief
 * @version 0.1
 * @date 2025-03-13
 *
 * @copyright
 *
 */

#ifndef PRUNING_HPP
#define PRUNING_HPP
#include<vector>
#include <grlina/graded_linalg.hpp>

namespace stable_decomposition {
    using index_t = int; // Change to large enough int type.

    using Mat = graded_linalg::R2GradedSparseMatrix<index_t>;
    /// Computes a pruning pair (I,K) of a module M given by presentation, following Bjerkevik 2025, Lemma 5.2
    std::vector<Mat> pruning_pair(const Mat &M, const double delta);
}

#endif // PRUNING_HPP