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

#include "include/pruning.hpp"
#include <thread>

namespace stable_decomposition {

using namespace graded_linalg;

template <typename Func, typename... Args>
auto timed_with_progress(const std::string &task_name, Func &&func,
                         Args &&...args)
    -> decltype(func(std::forward<Args>(args)...)) {

  std::atomic<bool> done(false);
  auto start = std::chrono::steady_clock::now();

  std::thread progress_thread([&]() {
    while (!done) {
      auto elapsed = std::chrono::steady_clock::now() - start;
      auto seconds =
          std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
      std::cout << "\r" << task_name << ": " << seconds << "s" << std::flush;
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
  });

  auto result = func(std::forward<Args>(args)...);

  done = true;
  progress_thread.join();

  auto end = std::chrono::steady_clock::now();
  auto total_seconds =
      std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
  std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Add this
  std::cout << "\r" << task_name << " completed in: " << total_seconds << "s"
            << std::endl;

  return result;
}

void matrix_reduction(vec<Mat>& A, vec<Mat>& B){
  index_t m = A.size();
  index_t n = B.size();
  std::unordered_map<index_t,index_t> pivots;
  vec<index_t> zero_columns = vec<index_t>();
  for (index_t j = 0; j < m+n; ++j) {
    index_t p;
    if(j < m){
      p = A[j].last_entry_index();
    } else {
      p = B[j-m].last_entry_index();
    }
    if (p != -1) {
      pivots[p] = j;
      for (index_t otherCol = 0; otherCol < m+n; ++otherCol) {
        if(otherCol < m){
          if (j != otherCol && A[otherCol].is_nonzero_at(p)) {
            if(j<m){
              add_to(A[j], A[otherCol]);
            } else {
              add_to( B[j-m], A[otherCol]);
            }
          }
        } else {
          if (j != otherCol && B[otherCol-m].is_nonzero_at(p)) {
            if(j<m){
              add_to(A[j], B[otherCol-m]);
            } else {
              add_to( B[j-m], B[otherCol-m]);
            }
          }
        }
      }
    } else {
      if(j >= m)
      zero_columns.push_back(j-m);
    }
  }
  vec_deletion(B, zero_columns); 
  // Todo: Probably should do this differently, since we expect there to be few non-zero matrices in B.
}

vec<Mat> homSpace(Mat& A, Mat& B) {
  // H: why does it allow me to compute rows forward when A is const where I'm
  // calling homSpace??
  // J: Your input is not Mat&A, but Mat A, so youre implicitely copying the input and then you can change it.
  // F: I don't understand the question, but the answer seems not to align with the function signature
  // J: seems to be fixed.
  A.compute_rows_forward();
  return hom_space_basis<r2degree, int, Mat>(
      A, B); // returns a basis of Hom(A, B) as a vector of matrices
}; // basis of Hom(A, B); done by hom_space_basis


vec<Mat> End_2d_0 (Mat &M, double delta) {
  Mat M_2d = M;
  M_2d.shift({delta, delta});//TODO F: Is this the shift by delta or 2*delta? // J: This is the shift by (delta, delta) in R^2
  vec<Mat> End_2d;
  try {
    End_2d = timed_with_progress("End_2d running", homSpace, M, M_2d);
    // std::cout << "Assignment completed" << std::endl; // Add this
  } catch (const std::exception &e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    throw;
  }
  auto dim_or = End_2d.size();
  vec<Mat> End_0;
  try {
    End_0 = timed_with_progress("End_0 running", homSpace, M, M);
    // std::cout << "Assignment completed" << std::endl; // Add this
  } catch (const std::exception &e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    throw;
  }

  matrix_reduction(End_0, End_2d);
  auto dim_quotient = End_2d.size();
  std::cout << "dim End_2d_0 = " << dim_quotient << " vs dim End_2d = " << dim_or << " and dim End_0 = " << End_0.size() << std::endl;
  return End_2d;
}

Mat zero_submodule(const Mat &m) {
  Mat zero(0, m.get_num_rows());      // Create a zero submodule of M
  zero.row_degrees = m.row_degrees;   // Ensure the row degrees match
  zero.col_degrees = vec<r2degree>(); // No columns, so no degrees
  zero.data = vec<vec<int>>();        // No data, so empty vector
  return zero;
}; // 0 as submodule of M, likely unnecessary

Mat all_submodule(const Mat &m) {
  Mat Id(m.get_num_rows(), m.get_num_rows(), "Identity");
  Id.row_degrees = m.row_degrees; // Ensure the row degrees match
  Id.col_degrees = m.row_degrees; // Ensure the column degrees match
  return Id;
}; // M as submodule of M, likely unnecessary

Mat submodule_sum(Mat l, Mat r) {
  assert(l.row_degrees == r.row_degrees);
  l.append_matrix(std::move(r));
  l.sort_columns_lexicographically();
  auto nzc = l.column_reduction_graded();
  l.delete_all_but_columns(nzc);
  return l;
}

//TODO H:  Check correctness.
Mat reduce_submodule(Mat M, Mat S) {//TODO Implementation looks like lots of copying, given that it's repeatedly called on the same matrix.
  auto M_copy = M;
  S.sort_columns_lexicographically();
  M_copy.append_matrix(S);
  auto nzc = M_copy.column_reduction_graded();
  int n = M.get_num_cols();
  for (int i = nzc.size() - 1; i >= 0; --i) {
    if (nzc[i] < n) {
      nzc.erase(nzc.begin() + i);
    }
  }
  if (nzc.size() == 0) {
    return zero_submodule(M);
  }
  M_copy.delete_all_but_columns(nzc);
  return M_copy;
}



Mat shifting_morphism(Mat A, double delta) {
  r2degree shift = {delta, delta};
  Mat result(A.get_num_rows(), A.get_num_rows(), "Identity");
  result.col_degrees = A.row_degrees;
  result.row_degrees = A.row_degrees;
  for (auto &g : result.row_degrees) {
    g = g - shift;
  }
  return result;
}

/// Return true of im M is contained in im N TODO: Is this correct? Relations
/// implicitely ignored.
bool image_contained_in_image(const Mat &M, const Mat &N) {
  auto N_copy = N;
  int threshold = N.get_num_cols();
  N_copy.append_matrix(M);
  N_copy.column_reduction_graded();
  bool is_zero = true;
  for (int i = threshold; i < N_copy.get_num_cols(); i++) {
    if (N_copy.data[i].size() != 0) {
      is_zero = false;
      break;
    }
  }
  return is_zero;
}

/// Checks if the two submodules generated by `A` and `B` of the module
/// presented by `M` are the same.
bool present_same_submodule(const Mat &M, const Mat &A, const Mat &B) {
  Mat MA = M, MB = M;
  MA.append_matrix(A);
  MB.append_matrix(B);
  return image_contained_in_image(B, MA) && image_contained_in_image(A, MB);
}

void print_progress(int iteration_I, size_t current, size_t total) {
  if(total != 0){
  int progress = (current * 100) / total;
  int bars = progress / 5; // 20 segments (100/5)

  std::cout << "\r Iteration " << iteration_I << " [";
  for (int i = 0; i < 20; ++i) {
    std::cout << (i < bars ? "█" : " ");
  }
  std::cout << "] " << progress << "%" << std::flush;
  }
}

std::pair<Mat, Mat> pruning_pair(Mat &M, const double delta) {
  // In this function, matrices either represent presentations of modules (relations -> generators),
  // or generators of a submodule (generators of submodule -> generators of module).
  // Therefore, matrices representing submodules only make sense w.r.t. an ambient module,
  // while a presentation of a module exists independently.

  Mat M2d = M;                                      // presentation matrix for M(2δ)
  M2d.shift({delta, delta});
  vec<Mat> B = End_2d_0(M, delta);                  // basis of Hom(M, M(2δ)).
  Mat can = shifting_morphism(M, delta);            // canonical morphism can: M → M(2δ)
  
  // Build the module I from the pruning pair (I,K)
  Mat I = all_submodule(M);                         // generators of Iᵢ ⊆ M, initialized I₀ = M ⊆ M
  int iteration_I = 1;
  for (;;) {
    assert(I.row_degrees == M.row_degrees);
    Mat I_new = I;                                  // Iᵢ₊₁ = ⋂_f f⁻¹(sh(Iᵢ, 2δ)),
    for (size_t idx = 0; const auto &f : B) {
      print_progress(iteration_I, ++idx, B.size());
      Mat foI = f * I_new;                          // generators of f(Iᵢ) ⊆ M(2δ)
      Mat canI = I_new;                             // generators of can(Iᵢ) ⊆ M(2δ)
                                                    //TODO F: Is it generators of can(I) \subseteq M or I \subseteq M(2d)?
      canI.shift_generators({delta, delta});
      auto inv = foI.inverse_image_copy(M2d, canI); // generators of f⁻¹(can(Iᵢ)) ⊆ I
      I_new = reduce_submodule(M, I_new * inv);     // generators of f⁻¹(can(Iᵢ)) ⊆ M
                                                    // H: without reduce_submodule, this gets suuuper slow
                                                    //TODO F: How can that be, shouldn't inverse image reduce?
    }
    print_progress(iteration_I, B.size(), B.size());
    std::cout << std::endl;
    if (present_same_submodule(M, I_new, I))        // run until stationary
      break;
    std::swap(I, I_new);
    iteration_I++;
  }

  // Build the module K from the pruning pair (I,K)
  auto canI = I;                                    // generators of can(I) ⊆ M(2δ)
  canI.shift_generators({delta, delta});
  Mat K = zero_submodule(M);                        // generators of Kᵢ ⊆ M, initialized K₀ = 0 ⊆ M
  int iteration_K = 1;
  for (;;) {
    assert(K.row_degrees == M.row_degrees);
    //TODO F: I changed the implementation slightly, because it was taking not sh⁻¹(f(Kᵢ), 2δ), but sh⁻¹(f([summation so far]), 2δ)
    //        I don't think that this has implications on correctness, but better reflects what we write I think.
    Mat K_new = zero_submodule(M);                  // Kᵢ₊₁ = ∑_f sh⁻¹(f(Kᵢ), 2δ)
    for (size_t idx = 0; const auto &f : B) {
      print_progress(iteration_K, ++idx, B.size());
      auto S = f * K_new;                           // generators of f(Kᵢ) ⊆ M(2δ)
      S = canI.inverse_image_copy(M2d, S);          // generators of sh⁻¹(f(Kᵢ), 2δ) ⊆ I
      S = I * S;                                    // generators of sh⁻¹(f(Kᵢ), 2δ) ⊆ M
      K_new.append_matrix(S);                       // K_{i+1} \coloneqq K_{}
    }
    K_new = reduce_submodule(M, K_new);
    print_progress(iteration_K, B.size(), B.size());
    std::cout << std::endl;
    if (present_same_submodule(M, K_new, K)) 
      break;
    std::swap(K, K_new);
    iteration_K++;
  }
  std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K" << std::endl;
  return {I, K};
}

Mat pruning(Mat &M, const double delta) {
  auto [I, K] = pruning_pair(M, delta);         // generators for I ⊆ M and K ⊆ M
  M.append_matrix(K);                           // presentation for M / K
  Mat Pru_M = I.presentation_of_submodule(M);   // presentation for I / K
  Pru_M.sort_columns_lexicographically();
  Pru_M.sort_rows_lexicographically();
  Pru_M.minimize();
  return Pru_M;
} 

std::optional<double> calculate_delta_from_matrix(const Mat& M) {
    const vec<r2degree>& col_degrees = M.col_degrees;
    const vec<r2degree>& row_degrees = M.row_degrees;
    
    if (col_degrees.empty() && row_degrees.empty()) {
        return std::nullopt;
    }
    
    // Initialize with first available degree
    double min_x, max_x, min_y, max_y;
    
    if (!col_degrees.empty()) {
        min_x = max_x = col_degrees[0].first;
        min_y = max_y = col_degrees[0].second;
    } else {
        min_x = max_x = row_degrees[0].first;
        min_y = max_y = row_degrees[0].second;
    }
    
    // Find bounding box of all degrees
    for (const auto& d : col_degrees) {
        min_x = std::min(min_x, d.first);
        max_x = std::max(max_x, d.first);
        min_y = std::min(min_y, d.second);
        max_y = std::max(max_y, d.second);
    }
    for (const auto& d : row_degrees) {
        min_x = std::min(min_x, d.first);
        max_x = std::max(max_x, d.first);
        min_y = std::min(min_y, d.second);
        max_y = std::max(max_y, d.second);
    }
    
    // Compute max extent
    double extent_x = max_x - min_x;
    double extent_y = max_y - min_y;
    double extent = std::max(extent_x, extent_y);
    
    if (extent < 1e-10) {
        return std::nullopt;
    }
    
    return extent * 0.01; // 1% of the extent
}

double get_delta(std::optional<double> user_delta, const Mat& M) {
    if (user_delta.has_value()) {
        std::cout << "Using specified delta: " << user_delta.value() << std::endl;
        return user_delta.value();
    }
    
    auto calculated = calculate_delta_from_matrix(M);
    if (calculated.has_value()) {
        std::cout << "Using calculated delta: " << calculated.value() << std::endl;
        return calculated.value();
    }
    
    std::cout << "Using default delta: 0.01" << std::endl;
    return 0.01;
}

}// namespace stable_decomposition

