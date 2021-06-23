cmake_minimum_required (VERSION 3.5...3.20)

project(PstClassifierSeqan CXX)

list(APPEND CMAKE_MODULE_PATH "/usr/share/cmake/Modules/") # For debian in singularity

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED CXX)
find_package(TBB REQUIRED)

find_package(Boost COMPONENTS thread system filesystem date_time iostreams REQUIRED)


set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/seqan3/build_system")
find_package(SeqAn3 REQUIRED)
set(SEQAN3_INCLUDED TRUE)

include_directories("${CMAKE_SOURCE_DIR}/robin-hood-hashing/src/include")

#set(Eigen3_DIR "${CMAKE_SOURCE_DIR}/eigen")
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories("${CMAKE_SOURCE_DIR}/eigen")

set(USE_BOOST OFF)
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN ON)
set(HIGHFIVE_EXAMPLES FALSE)

find_package(HDF5 REQUIRED COMPONENTS CXX)
add_subdirectory("${CMAKE_SOURCE_DIR}/HighFive")

add_subdirectory("${CMAKE_SOURCE_DIR}/kmc-vlmc")

# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

add_subdirectory(src)
add_subdirectory(tests)
