cmake_minimum_required(VERSION 3.30)
project(stable-decomposition)

set(CMAKE_CXX_STANDARD 23)
# Set details for release and debug mode. Use release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_DEBUG_POSTFIX _debug)

set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -flto") #TODO Want to compile with -Wall -Wextra -Wpedantic -Werror; currently too many warnings
set(CMAKE_CXX_FLAGS_DEBUG "-D_LIBCPP_ENABLE_ASSERTIONS=1 -D_GLIBCXX_DEBUG -g") # enable various sanity checks in the STL in libc++ and libstdc++, respectively.
set(SOURCES pruning.cpp)
find_package(Boost REQUIRED COMPONENTS timer graph )
include_directories(${Boost_INCLUDE_DIRS} Persistence-Algebra/include)

# Add executable defined in main.cpp
add_executable(pruning main.cpp ${SOURCES})
set_target_properties(pruning PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

# Add tests from tests/ subdirectory
include(CTest)
find_package(doctest REQUIRED)
include(doctest)
message(${CMAKE_CURRENT_LIST_DIR})
add_executable(tests tests/tests.cpp ${SOURCES})
target_link_libraries(tests doctest::doctest)
doctest_discover_tests(tests WORKING_DIRECTORY ../tests)