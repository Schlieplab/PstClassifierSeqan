cmake_minimum_required (VERSION 3.5...3.12)

project(PstClassifierSeqan CXX)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/seqan3/build_system")
find_package(SeqAn3 REQUIRED)
set(SEQAN3_INCLUDED TRUE)

include_directories("${CMAKE_SOURCE_DIR}/robin-hood-hashing/src/include")

#set(Eigen3_DIR "${CMAKE_SOURCE_DIR}/eigen")
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories("${CMAKE_SOURCE_DIR}/eigen")

set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN ON)
set(HIGHFIVE_EXAMPLES FALSE)
add_subdirectory("${CMAKE_SOURCE_DIR}/HighFive")


add_subdirectory(src)
add_subdirectory(tests)
