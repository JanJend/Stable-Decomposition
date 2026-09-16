#pragma once
#include <grlina/graded_linalg.hpp>

namespace stable_decomposition {
using index_t = int;
using Mat = graded_linalg::R2GradedSparseMatrix<index_t>;
using PModule = graded_linalg::R2Module<index_t>;
using OwnedSubmodule = graded_linalg::Submodule<Mat>;
using graded_linalg::vec;
} // namespace stable_decomposition
