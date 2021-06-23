cmake_minimum_required(VERSION 3.4...3.12)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts -fopenmp")

if (NOT SEQAN3_INCLUDED)
    set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
    find_package(SeqAn3 REQUIRED)
endif ()

find_package(Threads REQUIRED)

add_executable(pst-classifier pst-classifier.cpp)
target_link_libraries(pst-classifier seqan3::seqan3)
target_link_libraries(pst-classifier Threads::Threads)
target_link_libraries(pst-classifier CountVLMC::CountVLMC)

add_library(vlmc SHARED vlmc.cpp)
target_link_libraries(vlmc seqan3::seqan3)
target_link_libraries(vlmc Threads::Threads)
target_link_libraries(vlmc CountVLMC::CountVLMC)

add_executable(pst-batch-training batch-training.cpp)
target_link_libraries(pst-batch-training seqan3::seqan3)
target_link_libraries(pst-batch-training Threads::Threads)
target_link_libraries(pst-batch-training HighFive)
target_link_libraries(pst-batch-training ${HDF5_LIBRARIES})
target_link_libraries(pst-batch-training CountVLMC::CountVLMC)

add_executable(pst-score-sequences score-sequences.cpp)
target_link_libraries(pst-score-sequences seqan3::seqan3)
target_link_libraries(pst-score-sequences Threads::Threads)
target_link_libraries(pst-score-sequences HighFive)
target_link_libraries(pst-score-sequences ${HDF5_LIBRARIES})
target_link_libraries(pst-score-sequences CountVLMC::CountVLMC)


add_executable(bic bic.cpp)
target_link_libraries(bic seqan3::seqan3)
target_link_libraries(bic Threads::Threads)
target_link_libraries(bic CountVLMC::CountVLMC)

install(TARGETS pst-batch-training DESTINATION bin)
install(TARGETS pst-classifier DESTINATION bin)
install(TARGETS pst-score-sequences DESTINATION bin)

add_executable(composition-vectors composition-vectors.cpp)
target_link_libraries(composition-vectors seqan3::seqan3)
target_link_libraries(composition-vectors Threads::Threads)
target_link_libraries(composition-vectors HighFive)
target_link_libraries(composition-vectors ${HDF5_LIBRARIES})
target_link_libraries(composition-vectors CountVLMC::CountVLMC)

add_executable(calculate-distances calculate-distances.cpp)
target_link_libraries(calculate-distances seqan3::seqan3)
target_link_libraries(calculate-distances Threads::Threads)
target_link_libraries(calculate-distances HighFive)
target_link_libraries(calculate-distances ${HDF5_LIBRARIES})
target_link_libraries(calculate-distances CountVLMC::CountVLMC)

add_executable(frequency_differences frequency_differences.cpp)
target_link_libraries(frequency_differences seqan3::seqan3)
target_link_libraries(frequency_differences Threads::Threads)
target_link_libraries(frequency_differences ${HDF5_LIBRARIES})
target_link_libraries(frequency_differences CountVLMC::CountVLMC)
