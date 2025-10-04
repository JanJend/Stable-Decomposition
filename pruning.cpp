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

#include "pruning.hpp"
#include "grlina/sparse_matrix.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <thread>
#include <filesystem>
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
  // To-DO: Probably should do this differently, since we expect there to be few non-zero matrices in B.
}

vec<Mat> homSpace(Mat& A, Mat& B) {
  // H: why does it allow me to compute rows forward when A is const where I'm
  // calling homSpace??
  // J: Your input is not Mat&A, but Mat A, so youre implicitely copying the input and then you can change it.
  A.compute_rows_forward();
  return hom_space_basis<r2degree, int, Mat>(
      A, B); // returns a basis of Hom(A, B) as a vector of matrices
}; // basis of Hom(A, B); done by hom_space_basis


vec<Mat> End_2d_0 (Mat &M, double delta) {
  Mat M_2d = M;
  M_2d.shift(std::pair<double, double>(delta,delta));
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
  for(auto f : End_0){
    f.shift_generators(std::pair<double, double>(delta, delta));
  }
  matrix_reduction(End_0, End_2d);
  auto dim_quotient = End_2d.size();
  std::cout << "dim End_2d_0 = " << dim_quotient << " vs dim End_2d = " << dim_or << std::endl;
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
  auto nzc = l.column_reduction_graded();
  l.delete_all_but_columns(nzc);
  return l;
}

// H: Simplifies the presentation of S as a submodule of M. Check correctness.
Mat reduce_submodule(Mat M, Mat S) {
  auto M_copy = M;
  S.sort_columns_lexicographically();
  M_copy.append_matrix(S);
  M_copy.column_reduction_graded();
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



Mat image(const Mat &f, const Mat &A, const Mat &B,
          const Mat &U) { /// image of submodule U of A under morphism f: A -> B
  return f * U; // TODO not correct if submodules must be presented by minimal
                // generating systems
};


Mat shifting_morphism(Mat A, double delta) {
  r2degree shift = {delta, delta};
  Mat result(A.get_num_rows(), A.get_num_rows(), "Identity");
  result.col_degrees = A.row_degrees;
  result.row_degrees = A.row_degrees;
  for (auto &g : result.row_degrees) {
    g = g - shift;
  }
  return result;
  // return shifted_identity<r2degree, Mat>(A.row_degrees, shift);
}; // canonical morphism M -> M(2delta)

/// Return true of im M is contained in im N TO-DO: Is this correct? Relations
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
/// presented by `M` are isomorpic.
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



/// Computes a pruning pair (I,K) of a module M given by presentation, following
/// Bjerkevik 2025, Lemma 5.2
std::vector<Mat> pruning_pair(Mat &M, const double delta) {
  // Build a presentation matrix for M(2δ)
  Mat M2d = M;
  // H: Shift seems to work by shifting the grades up, which gives the opposite
  // of the positively shifted module
  // J: good catch, fixed that.
  M2d.shift({delta, delta});
  // Build a basis Γ for Hom(M, M(2δ)).
  vec<Mat> G = End_2d_0(M, delta);
  Mat can = shifting_morphism(M, delta);
  // Build the module I from the pruning pair
  Mat I = all_submodule(M);
  int iteration_I = 1;
  while (true) {
    Mat I_new = I;
    assert(I.get_num_rows() == M.get_num_rows());
    assert(I.row_degrees == M.row_degrees);
    // auto I_new = all_submodule(M);
    for (size_t idx = 0; idx < G.size(); ++idx) {
      print_progress(iteration_I, idx, G.size());
      const auto &f = G[idx];
      // (f \circ I_new)^{-1}(can*I)
      Mat foI = f * I_new;
      Mat I_shifted = I_new;
      // H: Using I_new instead of I_shifted = can * I_new; here is probably bad
      // form, since the row grades of I_new are wrong. But the row grades don't
      // matter, so this works (for now)
      // J: Fixed this to avoid unexpected behaviour. "can * -" should not be used, as it is inefficient.
      I_shifted.shift_generators(std::pair<double, double> (delta, delta));
      auto inv = foI.inverse_image_copy(M2d, I_shifted);
      //  H: without reduce_submodule, this gets suuuper slow
      I_new = reduce_submodule(M, I_new * inv);
      // J: Obsolete?
      // foI.append_matrix(I_new);
      // foI.append_matrix(M2d);
      // Mat K = foI.graded_kernel();
      // K.cull_columns(I_new.get_num_cols(), false);
      // K.column_reduction_graded();
      // I_new = M.submodule_intersection(I_new, f_copy.inverse_image(M2d,
      // I_shifted));
    }
    print_progress(iteration_I, G.size(), G.size());
    std::cout << std::endl;
    if (present_same_submodule(M, I_new, I)) {
      break;
    }
    std::swap(I, I_new);
    iteration_I++;
  }
  // Mat can_I; //TODO: 
  // the canonical morphism I -> I(2d)
  auto I_shifted = I;
  I_shifted.shift_generators(std::pair<double, double>(delta, delta));
  // Build the module K from the pruning pair
  Mat K = zero_submodule(M);
  int iteration_K = 1;
  while (true) {
    Mat K_new = K;
    // auto I_new = all_submodule(M);
    for (size_t idx = 0; idx < G.size(); ++idx) {
      assert(K.get_num_rows() == M.get_num_rows());
      print_progress(iteration_K, idx, G.size());
      const auto &f = G[idx];
      Mat K2 = K_new;
      auto f_K_new = f * K_new;
      K_new = I_shifted.inverse_image_copy(M2d, f_K_new);
      K_new = I * K_new;
      // This takes the sum with K from the previous step
      K_new.append_matrix(K2);
      K_new = reduce_submodule(M, K_new);
    }
    print_progress(iteration_K, G.size(), G.size());
    std::cout << std::endl;
    if (present_same_submodule(M, K_new, K)) {
      break;
    }
    std::swap(K, K_new);
    iteration_K++;
  }
  std::vector<Mat> pruning = {I, K};
  std::cout << iteration_I << " iterations for I, " << iteration_K << " iterations for K" << std::endl;
  
  
  return pruning; // Return the pruning pair (I, K)
}

} // namespace stable_decomposition

std::string insert_suffix_before_extension(const std::string& filepath, const std::string& suffix) {
    std::filesystem::path path(filepath);
    std::string stem = path.stem().string();             // filename without extension
    std::string extension = path.extension().string();   // e.g., ".txt"
    std::filesystem::path new_path = path.parent_path() / (stem + suffix + extension);
    return new_path.string();
}


int main(int argc, char** argv){
    std::string ex1 = "Persistence-Algebra/test_presentations/full_rips_size_1_instance_5_min_pres.scc";
    std::string ex2 = "tests/test5.scc";
    std::string torus = "tests/torus3_largestcomp.scc";
    std::string comp1 = "Persistence-Algebra/test_presentations/torus_100_0.10_dim1_decomposition/comp1.scc";
    std::string ex8 = "tests/test8.scc";
    std::string ex9 = "tests/test9.scc";
    std::string filepath;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file_path>" << std::endl;
        filepath = torus;
    } else {
        filepath = argv[1];
    }
    
    double delta = 0.01;
    std::filesystem::path input_path(filepath);
    std::cout << "Computing the pruning of " << input_path << std::endl;

    std::string modified_path = insert_suffix_before_extension(filepath, "_pru");
    std::filesystem::path output_path(modified_path);
    
    using namespace stable_decomposition;
    R2GradedSparseMatrix<int> M(filepath);
    
    M.sort_columns_lexicographically();
    M.sort_rows_lexicographically();
    // M.print_graded();
    std::vector<Mat> pruning = pruning_pair(M, delta);
    //print I
    //  std::cout << "Pruning pair (I, K):" << std::endl;
    // pruning[0].print_graded();
    //print K
    // pruning[1].print_graded();

    M.append_matrix(pruning[1]);
    Mat Pru_M = pruning[0].presentation_of_submodule(M);
    Pru_M.sort_columns_lexicographically();
    Pru_M.semi_minimize();
    std::ofstream output_file(output_path);
    if (!output_file.is_open()) {
        std::cerr << "Error: Unable to open output file " << output_path << std::endl;
        return 1;
    } else {
        Pru_M.to_stream(output_file);
        output_file.close();
        std::cout << "Pruning (" << delta << ") computed and saved to: " << output_path << std::endl;
    }
    
    return 0;
}
