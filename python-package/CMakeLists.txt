cmake_minimum_required (VERSION 3.4...3.12)
project(pst-classifier)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(Threads REQUIRED)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package(SeqAn3 REQUIRED)
set(SEQAN3_INCLUDED TRUE)

#set(Eigen3_DIR "${CMAKE_SOURCE_DIR}/../eigen/cmake")
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories("${CMAKE_SOURCE_DIR}/../robin-hood-hashing/src/include")
include_directories("${CMAKE_SOURCE_DIR}/../eigen")

set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN ON)
set(HIGHFIVE_EXAMPLES FALSE)
include_directories("${CMAKE_SOURCE_DIR}/../HighFive/include")
add_subdirectory("${CMAKE_SOURCE_DIR}/HighFive")

add_subdirectory(libvlmc)
