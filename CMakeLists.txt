cmake_minimum_required(VERSION 3.28)
project(stable-decomposition)
# C++ standard (overridable on the command line for testing newer standards).
# Kept at 17 to match AIDA, Persistence-Algebra and Skyscraper-Invariant.
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set details for release and debug mode. Use release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_DEBUG_POSTFIX _debug)
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -flto")
set(CMAKE_CXX_FLAGS_DEBUG "-D_LIBCPP_ENABLE_ASSERTIONS=1 -D_GLIBCXX_DEBUG -g")

# Source files
set(SOURCES 
    pruning.cpp
    utils.cpp
)

# Dependencies
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost COMPONENTS timer graph)
if(NOT ${Boost_FOUND})
    message(WARNING "Optional dependency Boost not found; please report if build errors occur.")
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
add_executable(pruning main.cpp ${SOURCES})
set_target_properties(pruning PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_link_libraries(pruning ${Boost_LIBRARIES})

# Tests (optional)
find_package(doctest QUIET)
if(doctest_FOUND)
    include(CTest)
    include(doctest)
    add_executable(tests tests/tests.cpp ${SOURCES})
    target_link_libraries(tests doctest::doctest ${Boost_LIBRARIES})
    doctest_discover_tests(tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
endif()