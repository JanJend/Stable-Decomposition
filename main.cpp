#include <string>
#include "pruning.hpp"
int main(){
    using namespace graded_linalg;
    using namespace stable_decomposition;
    std::string example_path1 = "Persistence-Algebra/test_presentations/full_rips_size_1_instance_5_min_pres.scc";
    std::string example_path2 = "test_presentations_pruning/test5.scc";
    R2GradedSparseMatrix<int> M(example_path2);
    M.print_graded();
    std::vector<Mat> pruning = pruning_pair(M, 1);
    //print I
    std::cout << "Pruning pair (I, K):" << std::endl;
    pruning[0].print_graded();
    //print K
    pruning[1].print_graded();
    return 0;
}