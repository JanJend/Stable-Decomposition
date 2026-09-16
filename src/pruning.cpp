

/**
 * @file pruning.cpp
 * @author Havard Bjerkevik, Jan Jendrysiak, and Fabian Lenzen
 * @brief
 * @version 0.2
 * @date 2025-10-6
 *
 * @copyright
 *
 */

#include "pruning.hpp"
#include <grlina/presentation_operations.hpp>
#include <numeric>
#include <algorithm>
#include <random>


namespace stable_decomposition {

using namespace graded_linalg;

std::pair<Mat, Mat> pruning_pair(Mat &M, const double delta, const bool quick) {
  // In this function, matrices either represent presentations of modules (relations -> generators),
  // or generators of a submodule (generators of submodule -> generators of module).
  // Therefore, matrices representing submodules only make sense w.r.t. an ambient module,
  // while a presentation of a module exists independently.
  // If quick==false, the number of iterations gives an upper bound on d_I(M,Pru_M).
  // If quick==true, we run a different version that often gives fewer iterations, but with no such bound

  Mat M2d = M;                                      // presentation matrix for M(2δ)
  M2d.shift({delta, delta});
  auto B = timed_with_progress("Shifted lifts", [&] {
    return graded_linalg::shifted_endomorphism_lift_complement(M, {delta, delta}, true);
  });

  /////////////////////////////////////////////////////////////////////
  ////////// Build the module I from the pruning pair (I,K) //////////
  /////////////////////////////////////////////////////////////////////
  // No additional maps impose constraints: I=M, K=0. In particular, the
  // quick-loop stopping criterion has no counter to visit when B is empty.
  if (B.empty()) return {graded_linalg::all_submodule(M), graded_linalg::zero_submodule(M)};
  Mat I = graded_linalg::all_submodule(M);                         // generators of Iᵢ ⊆ M, initialized I₀ = M ⊆ M
  int iteration_I = 1;
  std::vector<size_t> B_indices(B.size());
  B_indices.resize(B.size());                       // Gives the order we iterate over B. Randomized
  std::iota(B_indices.begin(), B_indices.end(), 0); // if quick==true
  if(quick){
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(B_indices.begin(), B_indices.end(), g);
    std::cout << std::endl;
  }
  size_t one_fifth = B.size();
  if(B.size() >= 20){
    one_fifth = B.size()/5;
  }
  Mat last_changed_I_new = I;
  size_t last_changed_iteration = B.size() - 1;
  bool done = false;
  Mat I_new = I;                                  // Iᵢ₊₁ = ⋂_f f⁻¹(sh(Iᵢ, 2δ)),
  for (;;) {
    assert(I.row_degrees == M.row_degrees);
    size_t idx = 0;
    Mat canI;
    if(!quick){
      canI = I_new;                               // generators of can(Iᵢ) ⊆ M(2δ)
                                                    //TODO F: Is it generators of can(I) \subseteq M or I \subseteq M(2d)?
      canI.shift_generators({delta, delta});
      canI = graded_linalg::reduce_submodule(M2d, canI);
    }
    //for (size_t idx = 0; const auto &f : B) {
    for (size_t counter = 0; counter < B.size(); counter++) {
      print_progress(iteration_I, ++idx, B.size());
      Mat& f = B[B_indices[counter]];
      Mat foI = f * I_new;                          // generators of f(Iᵢ) ⊆ M(2δ)
      if(quick){
        canI = I_new;                               // generators of can(Iᵢ) ⊆ M(2δ)
                                                      //TODO F: Is it generators of can(I) \subseteq M or I \subseteq M(2d)?
        canI.shift_generators({delta, delta});
        canI = graded_linalg::reduce_submodule(M2d, canI);
      }
      Mat inv = foI.inverse_image(M2d, canI); // generators of f⁻¹(can(Iᵢ)) ⊆ I
      Mat I_newxinv = I_new * inv;
      I_new = graded_linalg::reduce_submodule(M, I_newxinv);     // generators of f⁻¹(can(Iᵢ)) ⊆ M
                                                    // H: without reduce_submodule, this gets suuuper slow
                                                    //TODO F: How can that be, shouldn't inverse image reduce?
      // Break if we have gone one full run through all f without any change.
      // Check 5 times per run through B (5 is completely arbitrary)
      if(quick && (B.size() - 1 - counter) % one_fifth == 0){
        Mat MI_new = M;
        MI_new.append_matrix(I_new);
        if(!graded_linalg::image_contained_in_image(last_changed_I_new, MI_new)){
          last_changed_I_new = I_new;
          last_changed_iteration = counter;
          continue;
        }
        if(last_changed_iteration == counter){
          done = true;
          I = I_new;
          break;
        }
      }
    }

    print_progress(iteration_I, idx, B.size());
    std::cout << std::endl;
    if(quick){
      if(done){
        break;
      }
      iteration_I++;
      continue;
    }
    Mat MI_new = M;
    MI_new.append_matrix(I_new);
    if(graded_linalg::image_contained_in_image(I, MI_new)){
      break;
    }
    /* H: Could add the following check
    if (graded_linalg::image_contained_in_image(I_new, M)){
      Mat zro = graded_linalg::zero_submodule(M);
      return {zro, zro};
    }*/
    //if (graded_linalg::present_same_submodule(M, I_new, I))        // run until stationary
                                                    // H: I_new should automatically be a submod of I, so it
                                                    // should be enough to run image_contained_in_image once
                                                    // instead of twice, like it is in present_same_submodule
    //  break;
    I = I_new;
    iteration_I++;
  }

  /////////////////////////////////////////////////////////////////////
  ////////// Build the module K from the pruning pair (I,K) //////////
  /////////////////////////////////////////////////////////////////////
  done = false;
  Mat canI = I;                                    // generators of can(I) ⊆ M(2δ)
  canI.shift_generators({delta, delta});
  Mat K = graded_linalg::zero_submodule(M);                        // generators of Kᵢ ⊆ M, initialized K₀ = 0 ⊆ M
  Mat K_new = K;
  int iteration_K = 1;
  Mat last_changed_K_new = K;
  last_changed_iteration = B.size() - 1;
  for (;;) {
    size_t idx = 0;
    assert(K.row_degrees == M.row_degrees);
    //TODO F: I changed the implementation slightly, because it was taking not sh⁻¹(f(Kᵢ), 2δ), but sh⁻¹(f([summation so far]), 2δ)
    //        I don't think that this has implications on correctness, but better reflects what we write I think.
    if(!quick){
      Mat K_new = graded_linalg::zero_submodule(M);                  // Kᵢ₊₁ = ∑_f sh⁻¹(f(Kᵢ), 2δ)
    }
    for (size_t counter = 0; counter < B.size(); counter++) {
      print_progress(iteration_K, ++idx, B.size());
      Mat& f = B[B_indices[counter]];
      Mat S;
      if(quick){
        S = f * K_new;
      }else{
        S = f * K;                                    // generators of f(Kᵢ) ⊆ M(2δ)
      }
      canI = I;                                    // generators of can(I) ⊆ M(2δ)
      canI.shift_generators({delta, delta});
      S = canI.inverse_image(M2d, S);          // generators of sh⁻¹(f(Kᵢ), 2δ) ⊆ I
                                          // H: I'd like to use inverse_image_copy to avoid redefining
                                          // canI every loop, but iic doesn't compile.
      S = I * S;                                    // generators of sh⁻¹(f(Kᵢ), 2δ) ⊆ M
      K_new.append_matrix(S);                       // K_{i+1} \coloneqq K_{}
      K_new = graded_linalg::reduce_submodule(M, K_new);
      // Break if we have gone one full run through all f without any change.
      // Check 5 times per run through B (5 is completely arbitrary)
      if(quick && (B.size() - 1 - counter) % one_fifth == 0){
        Mat Mlast = M;
        Mlast.append_matrix(last_changed_K_new);
        if(!graded_linalg::image_contained_in_image(K_new, Mlast)){
          last_changed_K_new = K_new;
          last_changed_iteration = counter;
          continue;
        }
        if(last_changed_iteration == counter){
          done = true;
          K = K_new;
          break;
        }
      }
    }
    K_new = graded_linalg::reduce_submodule(M, K_new);
    print_progress(iteration_K, idx, B.size());
    std::cout << std::endl;
    if(quick){
      if(done){
        break;
      }
      iteration_K++;
      continue;
    }
    Mat MK = M;
    MK.append_matrix(K);
    if(graded_linalg::image_contained_in_image(K_new, MK))
    //if (graded_linalg::present_same_submodule(M, K_new, K))
      break;
    std::swap(K, K_new);
    iteration_K++;
  }
  std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K" << std::endl;
  return {I, K};
}

Mat pruning(Mat &M, const double delta, bool quick) {
  auto [I, K] = pruning_pair(M, delta, quick);         // generators for I ⊆ M and K ⊆ M
  //Mat K_module = K.presentation_of_submodule(M);
  //return K_module;
  M.append_matrix(K);                           // presentation for M / K
  Mat Pru_M = I.presentation_of_submodule(M);   // presentation for I / K
  Pru_M.sort_columns_lexicographically();
  Pru_M.sort_rows_lexicographically();
  Pru_M.column_reduction_graded_w_deletion();
  Pru_M.minimize();
  Pru_M.shift({-delta/2, -delta/2});
  return Pru_M;
}

std::pair<OwnedSubmodule, OwnedSubmodule> pruning_pair(
    std::shared_ptr<const PModule> module, const double delta, bool quick) {
  if (!module) throw std::invalid_argument("pruning_pair requires a module");
  Mat presentation = module->presentation();
  auto [I, K] = pruning_pair(presentation, delta, quick);
  return {OwnedSubmodule(module, std::move(I)),
          OwnedSubmodule(module, std::move(K))};
}

PModule pruning(PModule module, const double delta, bool quick) {
  Mat& presentation = module.mutable_presentation();
  return PModule(pruning(presentation, delta, quick));
}

} // namespace stable_decomposition
