cmake_minimum_required(VERSION 3.30)
project(stable-decomposition)

# Make produced makefile show full commands.
set(CMAKE_VERBOSE_MAKEFILE ON)
# Set details for release and debug mode. Use release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_DEBUG_POSTFIX _debug)
# Specify flags
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#TODO Want to compile with -Wall -Wextra -Wpedantic -Werror; currently to many warnings
set(CMAKE_CXX_FLAGS " -I${CMAKE_CURRENT_SOURCE_DIR}/Persistence-Algebra/include")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -flto")
# enable various sanity checks in the STL in libc++ and libstdc++, respectively. 
set(CMAKE_CXX_FLAGS_DEBUG "-D_LIBCPP_ENABLE_ASSERTIONS=1 -D_GLIBCXX_DEBUG -g")

set(SRC_FILES
	pruning.cpp
)

add_executable(pruning pruning.cpp ${SRC_FILES})
set_target_properties(pruning PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

find_package(Boost REQUIRED COMPONENTS timer graph )
include_directories(${Boost_INCLUDE_DIRS})

# find_package(OpenMP COMPONENTS CXX)
# if(OpenMP_CXX_FOUND)
#     message("Will compile with multiprocessing.")
#     target_compile_options(pruning PRIVATE ${OpenMP_CXX_FLAGS})
#     target_link_libraries(pruning ${OpenMP_CXX_FLAGS})
# else()
#     message("Will compile WITHOUT multiprocessing.")
# endif()

