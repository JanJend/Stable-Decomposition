

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
#include <grlina/hom_interface.hpp>
#include <numeric>
#include <algorithm>
#include <random>
#include <stdexcept>


namespace stable_decomposition {

using namespace graded_linalg;

// J: I am modularising the pruning function a bit for readability:

graded_linalg::r2degree pruning_shift(double epsilon) {
    if (!std::isfinite(2 * epsilon) || epsilon < 0 )
        throw std::invalid_argument("epsilon must be finite and nonnegative");
    return {2 * epsilon, 2 * epsilon};
}

/* New algorithm using the module framework*/
std::pair<Submodule, Submodule> pruning_pair(const Module& X, double epsilon, const bool quick, PruningProfile* profile) {

  const auto eps_vec = pruning_shift(epsilon);
  Hom eta_X = Homomorphism<Mat>::canonical_shift(X, eps_vec);

  auto B = measure(profile, &PruningProfile::basis, [&] {
    if (profile) return End_2d_0(eta_X, false);
    return timed_with_progress("End_2eps/0", [&] { return End_2d_0(eta_X, true); });
  });

 
  Submodule I = X.whole_submodule();
  // The following line is not consistent with the paper: 
  // I think that "K" should be Ker_2eps.
  if (B.empty()) return {I, X.zero_submodule()};

  // Compute I
  int iteration_I = 0;
  int counter = 0;
  {
  PruningTimer phase(profile ? &profile->i_phase : nullptr);
  Submodule I_ = X.zero_submodule();
  // J: The equality check takes some time here. Is it possible that in our example the submodules really dont change and we can just compare the matrices?
  while(!measure(profile, &PruningProfile::i_equal, [&] { return I.equals(I_); })){
    counter = 0;
    iteration_I++;
    if (profile) ++profile->i_iterations;
    I_ = I;
    auto J = measure(profile, &PruningProfile::i_image, [&] { return eta_X.image(I); });
    for(const Hom& f : B){
      if (!profile) print_progress(iteration_I, ++counter, B.size());
      auto preimage = measure(profile, &PruningProfile::i_preimage, [&] { return f.preimage(J); });
      I = measure(profile, &PruningProfile::i_intersection, [&] { return I.intersection(preimage, false); });
      measure(profile, &PruningProfile::i_reduce, [&] { I.reduce_generators_lazy(); });
    }
    
  }
  }
  if (!profile) print_progress(iteration_I, counter, B.size());

  // Compute K
  Submodule K = X.zero_submodule();
  int iteration_K = 0;
  
  {
  PruningTimer phase(profile ? &profile->k_phase : nullptr);
  Submodule K_ = X.whole_submodule();
  while(!measure(profile, &PruningProfile::k_equal, [&] { return K.equals(K_); })){
    counter = 0;
    iteration_K++;
    if (profile) ++profile->k_iterations;
    K_ = K;
    auto L = measure(profile, &PruningProfile::k_image, [&] { return eta_X.image(K); });
    for(const Hom& f : B){
      if (!profile) print_progress(iteration_K, ++counter, B.size());
      auto image = measure(profile, &PruningProfile::k_image, [&] { return f.image(K); });
      L = measure(profile, &PruningProfile::k_sum, [&] { return L.sum(image); });
    }
    // There is currently a difference here to the algorithm in the paper: 
    // eta_X i sused instead of eta_I, which means that K is not a submodule of I
    // I am mostly sure that this doesnt change anything.
    K = measure(profile, &PruningProfile::k_preimage, [&] { return eta_X.preimage(L); });
  }
  }
  if (!profile) print_progress(iteration_K, counter, B.size());

  if (!profile) std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K" << std::endl;
  if (!profile) std::cout << "Output is " << std::max(iteration_I, iteration_K)*epsilon << "-interleaved with the input" << std::endl;
  return {I,K};
}

Module pruning(Module X, const double epsilon, bool quick, PruningProfile* profile) {
  if (profile) *profile = {};
  PruningTimer total(profile ? &profile->total : nullptr);
  auto parent = std::make_shared<Module>(std::move(X));
  auto [I, K] = pruning_pair(*parent, epsilon, quick, profile);
  auto I_qK = measure(profile, &PruningProfile::quotient, [&I = I, &K = K] { return I.submodule_quotient(K); });
  measure(profile, &PruningProfile::presentation, [&] { I_qK.compute_presentation(); });
  Module result = std::move(I_qK);
  measure(profile, &PruningProfile::minimize, [&] { result.minimize(); });
  result.shift({-epsilon, -epsilon});
  return result;
}

std::pair<Mat, Mat> pruning_pair(Mat &M, const double epsilon, const bool quick, PruningProfile* profile) {
  // In this function, matrices either represent presentations of modules (relations -> generators),
  // or generators of a submodule (generators of submodule -> generators of module).
  // Therefore, matrices representing submodules only make sense w.r.t. an ambient module,
  // while a presentation of a module exists independently.
  // If quick==false, the number of iterations gives an upper bound on d_I(M,Pru_M).
  // If quick==true, we run a different version that often gives fewer iterations, but with no such bound

  const auto shift = pruning_shift(epsilon);
  Mat shifted_module = M;                                      // presentation matrix for M(2ε)
  shifted_module.shift(shift);
  auto B = measure(profile, &PruningProfile::basis, [&] {
    if (profile) return graded_linalg::End_2d_0(M, shift, false);
    return timed_with_progress("End_2eps/0", [&] { return graded_linalg::End_2d_0(M, shift, true); });
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
    if (!profile) std::cout << std::endl;
  }
  size_t one_fifth = B.size();
  if(B.size() >= 20){
    one_fifth = B.size()/5;
  }
  Mat last_changed_I_new = I;
  size_t last_changed_iteration = B.size() - 1;
  bool done = false;
  Mat I_new = I;                                  // Iᵢ₊₁ = ⋂_f f⁻¹(sh(Iᵢ, 2ε)),
  {
  PruningTimer phase(profile ? &profile->i_phase : nullptr);
  for (;;) {
    if (profile) ++profile->i_iterations;
    assert(I.row_degrees == M.row_degrees);
    size_t idx = 0;
    Mat canI;
    if(!quick){
      canI = I_new;                               // generators of can(Iᵢ) ⊆ M(2ε)
                                                    //TODO F: Is it generators of can(I) \subseteq M or I \subseteq M(2ε)?
      canI.shift_generators(shift);
      canI = measure(profile, &PruningProfile::i_reduce, [&] { return graded_linalg::reduce_submodule(shifted_module, canI); });
    }
    //for (size_t idx = 0; const auto &f : B) {
    for (size_t counter = 0; counter < B.size(); counter++) {
      if (!profile) print_progress(iteration_I, ++idx, B.size());
      Mat& f = B[B_indices[counter]];
      Mat foI = measure(profile, &PruningProfile::i_image, [&] { return f * I_new; });                          // generators of f(Iᵢ) ⊆ M(2ε)
      if(quick){
        canI = I_new;                               // generators of can(Iᵢ) ⊆ M(2ε)
                                                      //TODO F: Is it generators of can(I) \subseteq M or I \subseteq M(2ε)?
        canI.shift_generators(shift);
        canI = measure(profile, &PruningProfile::i_reduce, [&] { return graded_linalg::reduce_submodule(shifted_module, canI); });
      }
      Mat inv = measure(profile, &PruningProfile::i_preimage, [&] { return foI.inverse_image(shifted_module, canI); }); // generators of f⁻¹(can(Iᵢ)) ⊆ I
      Mat I_newxinv = measure(profile, &PruningProfile::i_image, [&] { return I_new * inv; });
      I_new = measure(profile, &PruningProfile::i_reduce, [&] { return graded_linalg::reduce_submodule(M, I_newxinv); });     // generators of f⁻¹(can(Iᵢ)) ⊆ M
                                                    // H: without reduce_submodule, this gets suuuper slow
                                                    //TODO F: How can that be, shouldn't inverse image reduce?
      // Break if we have gone one full run through all f without any change.
      // Check 5 times per run through B (5 is completely arbitrary)
      if(quick && (B.size() - 1 - counter) % one_fifth == 0){
        Mat MI_new = M;
        MI_new.append_matrix(I_new);
        if(!measure(profile, &PruningProfile::i_equal, [&] { return graded_linalg::image_contained_in_image(last_changed_I_new, MI_new); })){
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

    if (!profile) print_progress(iteration_I, idx, B.size());
    if (!profile) std::cout << std::endl;
    if(quick){
      if(done){
        break;
      }
      iteration_I++;
      continue;
    }
    Mat MI_new = M;
    MI_new.append_matrix(I_new);
    if(measure(profile, &PruningProfile::i_equal, [&] { return graded_linalg::image_contained_in_image(I, MI_new); })){
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

  }
  /////////////////////////////////////////////////////////////////////
  ////////// Build the module K from the pruning pair (I,K) //////////
  /////////////////////////////////////////////////////////////////////
  done = false;
  Mat canI = I;                                    // generators of can(I) ⊆ M(2ε)
  canI.shift_generators(shift);
  Mat K = graded_linalg::zero_submodule(M);                        // generators of Kᵢ ⊆ M, initialized K₀ = 0 ⊆ M
  Mat K_new = K;
  int iteration_K = 1;
  Mat last_changed_K_new = K;
  last_changed_iteration = B.size() - 1;
  {
  PruningTimer phase(profile ? &profile->k_phase : nullptr);
  for (;;) {
    if (profile) ++profile->k_iterations;
    size_t idx = 0;
    assert(K.row_degrees == M.row_degrees);
    //TODO F: I changed the implementation slightly, because it was taking not sh⁻¹(f(Kᵢ), 2ε), but sh⁻¹(f([summation so far]), 2ε)
    //        I don't think that this has implications on correctness, but better reflects what we write I think.
    if(!quick){
      Mat K_new = graded_linalg::zero_submodule(M);                  // Kᵢ₊₁ = ∑_f sh⁻¹(f(Kᵢ), 2ε)
    }
    for (size_t counter = 0; counter < B.size(); counter++) {
      if (!profile) print_progress(iteration_K, ++idx, B.size());
      Mat& f = B[B_indices[counter]];
      Mat S;
      if(quick){
        S = measure(profile, &PruningProfile::k_image, [&] { return f * K_new; });
      }else{
        S = measure(profile, &PruningProfile::k_image, [&] { return f * K; });                                    // generators of f(Kᵢ) ⊆ M(2ε)
      }
      canI = I;                                    // generators of can(I) ⊆ M(2ε)
      canI.shift_generators(shift);
      S = measure(profile, &PruningProfile::k_preimage, [&] { return canI.inverse_image(shifted_module, S); });          // generators of sh⁻¹(f(Kᵢ), 2ε) ⊆ I
                                          // H: I'd like to use inverse_image_copy to avoid redefining
                                          // canI every loop, but iic doesn't compile.
      S = measure(profile, &PruningProfile::k_image, [&] { return I * S; });                                    // generators of sh⁻¹(f(Kᵢ), 2ε) ⊆ M
      measure(profile, &PruningProfile::k_sum, [&] { K_new.append_matrix(S); });                       // K_{i+1} \coloneqq K_{}
      K_new = measure(profile, &PruningProfile::k_reduce, [&] { return graded_linalg::reduce_submodule(M, K_new); });
      // Break if we have gone one full run through all f without any change.
      // Check 5 times per run through B (5 is completely arbitrary)
      if(quick && (B.size() - 1 - counter) % one_fifth == 0){
        Mat Mlast = M;
        Mlast.append_matrix(last_changed_K_new);
        if(!measure(profile, &PruningProfile::k_equal, [&] { return graded_linalg::image_contained_in_image(K_new, Mlast); })){
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
    K_new = measure(profile, &PruningProfile::k_reduce, [&] { return graded_linalg::reduce_submodule(M, K_new); });
    if (!profile) print_progress(iteration_K, idx, B.size());
    if (!profile) std::cout << std::endl;
    if(quick){
      if(done){
        break;
      }
      iteration_K++;
      continue;
    }
    Mat MK = M;
    MK.append_matrix(K);
    if(measure(profile, &PruningProfile::k_equal, [&] { return graded_linalg::image_contained_in_image(K_new, MK); }))
    //if (graded_linalg::present_same_submodule(M, K_new, K))
      break;
    std::swap(K, K_new);
    iteration_K++;
  }
  }
  if (!profile) std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K" << std::endl;
  return {I, K};
}

Mat pruning(Mat &M, const double epsilon, bool quick, PruningProfile* profile) {
  if (profile) *profile = {};
  PruningTimer total(profile ? &profile->total : nullptr);
  auto [I, K] = pruning_pair(M, epsilon, quick, profile);         // generators for I ⊆ M and K ⊆ M
  //Mat K_module = K.presentation_of_submodule(M);
  //return K_module;
  measure(profile, &PruningProfile::quotient, [&M, &K = K] { M.append_matrix(K); });                           // presentation for M / K
  Mat Pru_M = measure(profile, &PruningProfile::presentation, [&I = I, &M] { return I.presentation_of_submodule(M); });   // presentation for I / K
  measure(profile, &PruningProfile::minimize, [&] {
    Pru_M.sort_columns_lexicographically();
    Pru_M.sort_rows_lexicographically();
    Pru_M.column_reduction_graded_w_deletion();
    Pru_M.minimize();
  });
  Pru_M.shift({-epsilon, -epsilon});
  return Pru_M;
}






} // namespace stable_decomposition
