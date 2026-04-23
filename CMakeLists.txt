cmake_minimum_required(VERSION 3.28)
project(stable-decomposition)
set(CMAKE_CXX_STANDARD 23)

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
cmake_policy(SET CMP0167 NEW) # Use new mechanism for finding Boost
find_package(Boost REQUIRED COMPONENTS timer graph)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} Persistence-Algebra/include) # suppress warnings for system includes

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