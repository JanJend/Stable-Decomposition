
#include "include/grlina/r2graded_matrix.hpp"  
#include "include/grlina/r3graded_matrix.hpp"  
#include <iostream>
#include <cassert>

using namespace graded_linalg; 

void test() {
    R2GradedSparseMatrix<int> M("/home/wsljan/AIDA/test_presentations/full_rips_size_1_instance_5_min_pres.scc");
    M.print_graded();
    M.sort_columns_lexicographically();
    M.print_graded();
}

int main() {
    test();
    return 0;
}