cmake_minimum_required(VERSION 3.15)
project(pstclassifierseqan)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(src/search)
include_directories(src/search/lazy_suffix_tree)
include_directories(tests)

add_executable(pstclassifierseqan
        src/search/lazy_suffix_tree/construction.hpp
        src/search/lazy_suffix_tree/iteration.hpp
        src/search/lazy_suffix_tree/suffix_links.hpp
        src/search/lazy_suffix_tree.hpp
        src/CMakeLists.txt
        src/probabilistic_suffix_tree.hpp
        src/pst-classifier.cpp
        src/vlmc.cpp
        tests/search/lazy_suffix_tree/construction_tests.cpp
        tests/search/lazy_suffix_tree/iteration_tests.cpp
        tests/search/lazy_suffix_tree/suffix_links_tests.cpp
        tests/search/lazy_suffix_tree_tests.cpp
        tests/benchmarks.cpp
        tests/probabilistic_suffix_tree_tests.cpp
        CMakeLists.txt)
