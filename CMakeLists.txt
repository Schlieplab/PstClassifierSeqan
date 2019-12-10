cmake_minimum_required (VERSION 3.4...3.12)

project(PstClassifierSeqan CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/seqan3/build_system")
find_package(SeqAn3 REQUIRED)
set(SEQAN3_INCLUDED TRUE)

add_subdirectory(src)
add_subdirectory(tests)
