#include "algebra_compat.hpp"
#include "progress.hpp"
#include <grlina/module/submodule.hpp>
#include <grlina/hom_interface.hpp>
#include <grlina/matrices/matrix_algorithms.hpp>

namespace stable_decomposition {
void matrix_reduction(vec<Mat>& A, vec<Mat>& B) {
    graded_linalg::reduce_matrix_family_modulo(A, B);
}
vec<Mat> homSpace(Mat& A, Mat& B) {
    // Preserve the legacy row-cache side effect for callers using it afterward.
    A.compute_rows_forward();
    return graded_linalg::homomorphism_lift_basis(A, B);
}
vec<Mat> End_2epsilon_0(Mat& M, double epsilon) {
    return timed_with_progress("End_2eps/0", [&] {
        return graded_linalg::End_2d_0(M, pruning_shift(epsilon), true);
    });
}
vec<Mat> End_2d_0(Mat& M, double epsilon) { return End_2epsilon_0(M, epsilon); }
Mat zero_submodule(const Mat& M) { return graded_linalg::zero_submodule(M); }
Mat all_submodule(const Mat& M) { return graded_linalg::all_submodule(M); }
Mat submodule_sum(Mat A, Mat B) { return graded_linalg::submodule_sum(A, B); }
Mat reduce_submodule(Mat& M, Mat& S) { return graded_linalg::reduce_submodule(M, S); }
Mat shifting_morphism(Mat A, double amount) {
    return graded_linalg::canonical_shift_lift(A, {amount, amount});
}
bool image_contained_in_image(const Mat& A, const Mat& B) {
    Mat presentation = graded_linalg::empty_matrix_like(B, 0, B.get_num_rows());
    presentation.row_degrees = B.row_degrees;
    auto parent = std::make_shared<const Module>(std::move(presentation));
    return Submodule(parent, B).contains(Submodule(parent, A));
}
bool present_same_submodule(const Mat& M, const Mat& A, const Mat& B) {
    return graded_linalg::present_same_submodule(M, A, B);
}
Mat image(const Mat& f, const Mat& A, const Mat& B, const Mat& U) {
    auto domain = std::make_shared<const Module>(A);
    auto target = std::make_shared<const Module>(B);
    return graded_linalg::Homomorphism<Mat>(domain, target, f)
        .image(Submodule(domain, U)).generator_map().generator_lift();
}
} // namespace stable_decomposition
