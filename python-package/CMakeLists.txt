cmake_minimum_required (VERSION 3.4...3.12)
project(pst-classifier)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(Threads REQUIRED)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package(SeqAn3 REQUIRED)
set(SEQAN3_INCLUDED TRUE)


include_directories("${CMAKE_SOURCE_DIR}/../robin-hood-hashing/src/include")
include_directories("${CMAKE_SOURCE_DIR}/../eigen")

add_subdirectory(libvlmc)
