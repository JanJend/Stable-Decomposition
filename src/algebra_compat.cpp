#include "algebra_compat.hpp"
#include "progress.hpp"
#include <grlina/presentation_operations.hpp>

namespace stable_decomposition {
void matrix_reduction(vec<Mat>& A, vec<Mat>& B) {
    graded_linalg::reduce_matrix_family_modulo(A, B);
}
vec<Mat> homSpace(Mat& A, Mat& B) {
    // Preserve the legacy row-cache side effect for callers using it afterward.
    A.compute_rows_forward();
    return graded_linalg::homomorphism_lift_basis(A, B);
}
vec<Mat> End_2d_0(Mat& M, double delta) {
    return timed_with_progress("Shifted lifts", [&] {
        return graded_linalg::shifted_endomorphism_lift_complement(M, {delta, delta}, true);
    });
}
Mat zero_submodule(const Mat& M) { return graded_linalg::zero_submodule(M); }
Mat all_submodule(const Mat& M) { return graded_linalg::all_submodule(M); }
Mat submodule_sum(Mat A, Mat B) { return graded_linalg::submodule_sum(A, B); }
Mat reduce_submodule(Mat& M, Mat& S) { return graded_linalg::reduce_submodule(M, S); }
Mat shifting_morphism(Mat A, double delta) {
    return graded_linalg::canonical_shift_lift(A, {delta, delta});
}
bool image_contained_in_image(const Mat& A, const Mat& B) {
    return graded_linalg::image_contained_in_image(A, B);
}
bool present_same_submodule(const Mat& M, const Mat& A, const Mat& B) {
    return graded_linalg::present_same_submodule(M, A, B);
}
Mat image(const Mat& f, const Mat& A, const Mat& B, const Mat& U) {
    auto domain = std::make_shared<const PModule>(A);
    auto target = std::make_shared<const PModule>(B);
    return graded_linalg::Homomorphism<Mat>(domain, target, f)
        .image(OwnedSubmodule(domain, U)).generators();
}
} // namespace stable_decomposition
