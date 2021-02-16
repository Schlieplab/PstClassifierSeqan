cmake_minimum_required (VERSION 3.4...3.12)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")


if (NOT SEQAN3_INCLUDED)
    set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
    find_package (SeqAn3 REQUIRED)
endif()

find_package(Threads REQUIRED)

add_executable (pst-classifier pst-classifier.cpp)
target_link_libraries (pst-classifier seqan3::seqan3)
target_link_libraries(pst-classifier Threads::Threads)

add_library (vlmc SHARED vlmc.cpp)
target_link_libraries (vlmc seqan3::seqan3)
target_link_libraries(vlmc Threads::Threads)

add_executable (batch-training batch-training.cpp)
target_link_libraries (batch-training seqan3::seqan3)
target_link_libraries(batch-training Threads::Threads)
target_link_libraries(batch-training HighFive)

add_executable (score-sequences score-sequences.cpp)
target_link_libraries (score-sequences seqan3::seqan3)
target_link_libraries(score-sequences Threads::Threads)
target_link_libraries(score-sequences HighFive)