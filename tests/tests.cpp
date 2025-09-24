#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <string>
#include "../pruning.hpp"
bool run_tests(std::string infile){
    using namespace stable_decomposition;
    Mat M(infile);
    std::vector<Mat> pruning = pruning_pair(M, 1);
    return true;
}
TEST_CASE("test1.scc") { run_tests("test1.scc"); }
TEST_CASE("test2.scc") { run_tests("test2.scc"); }
TEST_CASE("test3.scc") { run_tests("test3.scc"); }
TEST_CASE("test4.scc") { run_tests("test4.scc"); }
TEST_CASE("test5.scc") { run_tests("test5.scc"); }
