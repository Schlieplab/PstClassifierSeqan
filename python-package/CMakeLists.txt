cmake_minimum_required(VERSION 3.11.0)
project(pst-classifier)


find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)


set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package(SeqAn3 REQUIRED)
set(SEQAN3_INCLUDED TRUE)


add_subdirectory(libvlmc)
