cmake_minimum_required (VERSION 3.4...3.12)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if (NOT SEQAN3_INCLUDED)
    set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
    find_package (SeqAn3 REQUIRED)
endif()

add_executable (pst-classifier pst-classifier.cpp)
target_link_libraries (pst-classifier seqan3::seqan3)

add_library (vlmc SHARED vlmc.cpp)
target_link_libraries (vlmc seqan3::seqan3)
