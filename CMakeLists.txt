cmake_minimum_required(VERSION 3.28)
project(stable-decomposition VERSION 0.3.0 LANGUAGES CXX)
# C++ standard (overridable on the command line for testing newer standards).
# Kept at 17 to match AIDA and Persistence-Algebra.
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Keep storage consistent in the CLI, tests, and embedded AIDA sources.
option(GRLINA_USE_COMPACT_COLUMNS "Use compact columns instead of std::vector columns" ON)
add_compile_definitions(GRLINA_USE_COMPACT_COLUMNS=$<BOOL:${GRLINA_USE_COMPACT_COLUMNS}>)

# Set details for release and debug mode. Use release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_DEBUG_POSTFIX _debug)
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -flto")
if(APPLE)
    set(CMAKE_CXX_FLAGS_DEBUG "-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG -g")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-D_GLIBCXX_DEBUG -g")
endif()

# Source files
set(SOURCES 
    src/pruning.cpp
    src/pruning_profiled.cpp
    src/utils.cpp
    src/algebra_compat.cpp
    src/epsilon.cpp
    src/progress.cpp
)

# Dependencies
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED COMPONENTS timer graph)

# Reuse AIDA's parser dependency and call its library directly when available.
set(AIDA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../AIDA" CACHE PATH "AIDA checkout")
find_path(CLI11_INCLUDE_DIR CLI11.hpp HINTS "${AIDA_DIR}/include" REQUIRED)
option(PRUNING_WITH_AIDA "Enable direct AIDA decomposition comparison" ON)
if(PRUNING_WITH_AIDA)
    add_subdirectory("${AIDA_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/aida" EXCLUDE_FROM_ALL)
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
# Persistence-Algebra is resolved as a sibling checkout, the same way AIDA does it,
# so all repos share one working copy instead of a separate nested submodule.
set(PERSISTENCE_ALGEBRA_DIR "${CMAKE_SOURCE_DIR}/../Persistence-Algebra" CACHE PATH "Path to Persistence-Algebra directory")
if(NOT EXISTS "${PERSISTENCE_ALGEBRA_DIR}/include/grlina/graded_linalg.hpp")
    message(FATAL_ERROR "Persistence-Algebra not found at ${PERSISTENCE_ALGEBRA_DIR}. "
                        "Clone it next to this repository or set PERSISTENCE_ALGEBRA_DIR.")
endif()
# suppress warnings for system includes
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${PERSISTENCE_ALGEBRA_DIR}/include)

# OpenMP is optional: used when available, ignored otherwise.
# Apple clang ships no OpenMP runtime, so point CMake at Homebrew's libomp if present.
if(APPLE AND NOT DEFINED OpenMP_ROOT)
    execute_process(COMMAND brew --prefix libomp
                    OUTPUT_VARIABLE _brew_libomp
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    if(_brew_libomp AND EXISTS "${_brew_libomp}")
        set(OpenMP_ROOT "${_brew_libomp}")
    endif()
endif()
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
endif()

# Main executable
add_library(stable_pruning STATIC ${SOURCES})
find_package(Threads REQUIRED)
target_link_libraries(stable_pruning PUBLIC ${Boost_LIBRARIES} Threads::Threads)
target_include_directories(stable_pruning SYSTEM PRIVATE "${CLI11_INCLUDE_DIR}")
target_compile_definitions(stable_pruning PRIVATE PRUNING_VERSION="${PROJECT_VERSION}")
target_compile_definitions(stable_pruning PUBLIC GRLINA_USE_COMPACT_COLUMNS=$<BOOL:${GRLINA_USE_COMPACT_COLUMNS}>)
add_library(pruning_cli STATIC src/cli.cpp)
target_link_libraries(pruning_cli PUBLIC stable_pruning)
if(APPLE)
    target_compile_definitions(pruning_cli PRIVATE GRLINA_HILBERT_CORETEXT)
    target_link_libraries(pruning_cli PRIVATE "-framework CoreText" "-framework CoreGraphics" "-framework CoreFoundation")
endif()
if(PRUNING_WITH_AIDA)
    target_link_libraries(pruning_cli PRIVATE aida_static)
    target_include_directories(pruning_cli BEFORE PRIVATE "${AIDA_DIR}/include" "${AIDA_DIR}/src")
    target_compile_definitions(pruning_cli PUBLIC PRUNING_WITH_AIDA=1)

    # Explicit target: experiments are never part of the default build or CTest.
    add_executable(pruning_experiment EXCLUDE_FROM_ALL experiments/support/runtime.cpp)
    target_link_libraries(pruning_experiment PRIVATE stable_pruning aida_static)
    target_include_directories(pruning_experiment BEFORE PRIVATE "${AIDA_DIR}/include" "${AIDA_DIR}/src")
endif()
add_executable(pruning main.cpp)
set_target_properties(pruning PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_link_libraries(pruning PRIVATE pruning_cli)

# Explicit benchmark target; no AIDA dependency or generated output modules.
add_executable(pruning_optimization_worker EXCLUDE_FROM_ALL experiments/support/optimization_runtime.cpp)
target_link_libraries(pruning_optimization_worker PRIVATE stable_pruning)

# Tests (optional)
include(CTest)
if(BUILD_TESTING)
    add_executable(pruning_optimizations_test tests/functional/pruning_optimizations_test.cpp)
    target_link_libraries(pruning_optimizations_test PRIVATE stable_pruning)
    target_compile_definitions(pruning_optimizations_test PRIVATE
        PRUNING_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_presentations"
        PRUNING_TEST_PRESENTATIONS_DIR="${PERSISTENCE_ALGEBRA_DIR}/test_presentations")
    add_test(NAME pruning_optimizations COMMAND pruning_optimizations_test)
    set_tests_properties(pruning_optimizations PROPERTIES TIMEOUT 120)
    add_executable(pruning_solutions_test tests/pruning_solutions_test.cpp)
    target_link_libraries(pruning_solutions_test PRIVATE stable_pruning)
    target_compile_definitions(pruning_solutions_test PRIVATE
        PRUNING_SOLUTIONS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_presentations/tests_with_solutions")
    add_test(NAME pruning_solutions COMMAND pruning_solutions_test)
    set_tests_properties(pruning_solutions PROPERTIES TIMEOUT 120)
    add_executable(csc_pruning_boundary_test tests/functional/csc_pruning_boundary_test.cpp)
    target_link_libraries(csc_pruning_boundary_test PRIVATE stable_pruning)
    target_compile_options(csc_pruning_boundary_test PRIVATE -UNDEBUG)
    add_test(NAME csc_pruning_boundary COMMAND csc_pruning_boundary_test)
    add_executable(hilbert_plot_test tests/functional/hilbert_plot_test.cpp)
    add_test(NAME hilbert_plot COMMAND hilbert_plot_test)
    add_test(NAME pruning_cli_outputs COMMAND ${CMAKE_COMMAND}
        -DPRUNING=$<TARGET_FILE:pruning> -DWITH_AIDA=${PRUNING_WITH_AIDA} -DVERSION=${PROJECT_VERSION}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/functional/cli_outputs.cmake)
    set_tests_properties(pruning_cli_outputs PROPERTIES TIMEOUT 180)
    add_executable(compare_pruning tests/functional/compare_pruning.cpp)
    target_link_libraries(compare_pruning PRIVATE stable_pruning)
    target_compile_definitions(compare_pruning PRIVATE
        PRUNING_TEST_PRESENTATIONS_DIR="${PERSISTENCE_ALGEBRA_DIR}/test_presentations")
    if(EXISTS "${PERSISTENCE_ALGEBRA_DIR}/test_presentations/toy_example_1.scc")
        add_test(NAME compare_pruning_zero COMMAND compare_pruning --epsilon 0 toy_example_2.scc)
        add_test(NAME compare_pruning_examples COMMAND compare_pruning)
        set_tests_properties(compare_pruning_zero compare_pruning_examples PROPERTIES TIMEOUT 120)
    endif()
    if(EXISTS "${PERSISTENCE_ALGEBRA_DIR}/test_presentations/small_pres_with_nonzero_k_min_pres.firep")
        add_test(NAME compare_pruning_nonzero_k COMMAND compare_pruning --epsilon 0.5
            small_pres_with_nonzero_k_min_pres.firep)
        set_tests_properties(compare_pruning_nonzero_k PROPERTIES TIMEOUT 120)
    endif()
endif()
find_package(doctest QUIET)
if(BUILD_TESTING AND doctest_FOUND)
    include(CTest)
    include(doctest)
    add_executable(tests tests/functional/tests.cpp tests/functional/algebra_adapters_test.cpp tests/functional/epsilon_test.cpp tests/functional/dynamic_grid_pruning_test.cpp)
    target_link_libraries(tests PRIVATE doctest::doctest stable_pruning)
    target_compile_definitions(tests PRIVATE
        PRUNING_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_presentations")
    doctest_discover_tests(tests)
    add_test(NAME pruning_requires_input COMMAND $<TARGET_FILE:pruning>)
    set_tests_properties(pruning_requires_input PROPERTIES WILL_FAIL TRUE)
    add_test(NAME pruning_epsilon_cli COMMAND ${CMAKE_COMMAND}
        -DPRUNING=$<TARGET_FILE:pruning> -DINPUT=${CMAKE_SOURCE_DIR}/test_presentations/test1.scc
        -P ${CMAKE_SOURCE_DIR}/tests/functional/epsilon_cli.cmake)
    set_tests_properties(pruning_epsilon_cli PROPERTIES TIMEOUT 120)
endif()
