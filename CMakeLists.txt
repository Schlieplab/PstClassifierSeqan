cmake_minimum_required (VERSION 3.4)
project(PstClassifierSeqan)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)
set(SEQAN3_INCLUDED TRUE)
add_subdirectory(tests)
