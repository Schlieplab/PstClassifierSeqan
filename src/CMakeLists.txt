cmake_minimum_required (VERSION 3.4)
project (pst_classifier CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package (SeqAn3 REQUIRED)

add_executable (vlmc vlmc.cpp)

target_link_libraries (vlmc seqan3::seqan3)
