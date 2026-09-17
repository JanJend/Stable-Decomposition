#pragma once
#include <grlina/graded_linalg.hpp>

namespace stable_decomposition {
using index_t = int;
using Mat = graded_linalg::R2GradedSparseMatrix<index_t>;
using Module = graded_linalg::R2Module<index_t>;
using Submodule = graded_linalg::Submodule<Mat>;
using graded_linalg::vec;
using Hom = graded_linalg::Homomorphism<Mat>;

} // namespace stable_decomposition
