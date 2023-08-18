cmake_minimum_required (VERSION 3.4...3.12)
project(pst-classifier)

list(APPEND CMAKE_MODULE_PATH "/usr/share/cmake/Modules/") # For debian in singularity
list(APPEND CMAKE_MODULE_PATH "/usr/lib/python3.11/site-packages/skbuild/resources/cmake/") # For skbuild

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package(SeqAn3 REQUIRED)
set(SEQAN3_INCLUDED TRUE)


set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN ON)
set(HIGHFIVE_EXAMPLES FALSE)

include_directories("${CMAKE_SOURCE_DIR}/../HighFive/include")

include_directories("${CMAKE_SOURCE_DIR}/../robin-hood-hashing/src/include")
include_directories("${CMAKE_SOURCE_DIR}/../eigen")
include_directories("${CMAKE_SOURCE_DIR}/../indicators/include")

add_subdirectory("${CMAKE_SOURCE_DIR}/../vlmc-from-kmers" "vlmc-from-kmers")
include_directories("${CountVLMC_INCLUDE_DIRS}")

include_directories("${CMAKE_SOURCE_DIR}/../include/")

add_subdirectory(libvlmc)
