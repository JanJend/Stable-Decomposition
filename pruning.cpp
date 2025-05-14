/**
 * @file pruning.cpp
 * @author Jan Jendrysiak and Fabian Lenzen
 * @brief
 * @version 0.1
 * @date 2025-03-13
 *
 * @copyright
 *
 */

#include "grlina/graded_matrix.hpp"
#include <grlina/graded_linalg.hpp>
#include <algorithm>
namespace stable_decomposition {

using namespace graded_linalg;
using degree_t = r2degree;
using index_t = int; // Change to large enough int type.
using Mat = GradedSparseMatrix<degree_t, index_t>;

Mat submodule_intersection(Mat l, Mat r){
    assert(l.row_degrees == r.row_degrees);
    auto lr = l;
    lr.append_matrix(std::move(r));
    auto k = lr.graded_kernel();
    k.cull_columns(r.get_num_cols());
    return l*k;//TODO implement
};

Mat submodule_sum(Mat l, Mat r){
    assert(l.row_degrees == r.row_degrees);
    l.append_matrix(std::move(r));
    auto nzc = l.column_reduction_graded();
    l.delete_all_but_columns(nzc);
    return l;
}

Mat image(Mat f, Mat U){// image of U \subseteq A under morphism f: A -> B

};
Mat preimage(Mat f, Mat U);// preimage of U \subseteq B under morphism f: A -> B
vec<Mat> hom_space(Mat A, Mat B); // basis of Hom(A, B);
Mat shifting_morphism(Mat A, double delta); // canonical morphism M -> M(2delta)
Mat zero_submodule(Mat m); // 0 as submodule of M, likely unnecessary
Mat all_submodule(Mat m); // M as submodule of M, likely unnecessary

/// Computes a pruning pair (I,K) of a module M given by presentation, following Bjerkevik 2025, Lemma 5.2
auto pruning_pair(GradedSparseMatrix<degree_t, index_t> M, const double delta){
    // Build a presentation matrix for N := M(2δ)
    auto N = M;
    for(auto &&g : N.row_degrees) {g.first += delta; g.second += delta;}
    for(auto &&g : N.col_degrees) {g.first += delta; g.second += delta;}
    // Build a basis Γ for Hom(M, M(2δ)).
    vec<SparseMatrix<index_t>> G; //TODO: make G a basis of hom_space(M, N);
    GradedSparseMatrix<degree_t, index_t> can; //TODO: build the canonical morphism M -> M(2d)
    // Build the module I from the pruning pair
    auto I = M, Ii = M;
    for(;;) {
        auto Ii_new = M;
        for(auto &&f : G){
            // TODO need module intersection, preimage and image under module homomorphism
            Ii_new = Ii_new.intersect(preimage(f, image(can, Ii)));
        }
        auto I_new = I.intersect(Ii_new);
        if(presents_same_module(I_new, I)) //TODO implement
            break;
        I = std::move(I_new);
        Ii = std::move(Ii_new);
    }
    GradedSparseMatrix<degree_t, index_t> can_I; //TODO: the canonical morphism I -> I(2d)
    // Build the module K from the pruning pair
    auto K = zero_submodule(M), Ki = can_I.graded_kernel(); //TODO implement
    for(;;){
        auto Ki_new = zero_submodule(M);
        for(auto &&f : G){
            Ki_new = Ki_new.submodule_sum(preimage(can_I, image(f, Ki)));
        }
        auto K_new = K.submodule_sum(Ki_new);
        if(presents_same_module(K_new, K)) //TODO implement
            break;
        K = std::move(K_new);
        Ki = std::move(Ki_new);
    }

}
} // namespace stable_decomposition
