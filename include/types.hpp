#pragma once
#include <grlina/graded_linalg.hpp>

namespace stable_decomposition {
using index_t = int;
using Mat = graded_linalg::R2GradedSparseMatrix<index_t>;
inline constexpr const char* column_storage =
    graded_linalg::is_compact_storage_v<decltype(Mat::data)> ? "compact" : "vector";
using Module = graded_linalg::Module<Mat>;
using Submodule = graded_linalg::Submodule<Mat>;
using graded_linalg::vec;
using Hom = graded_linalg::Homomorphism<Mat>;

} // namespace stable_decomposition
