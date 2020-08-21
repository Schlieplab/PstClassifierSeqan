cmake_minimum_required (VERSION 3.4...3.12)

project (pst_classifier_tests CXX)


include(CTest)
find_package(Threads REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

if (NOT SEQAN3_INCLUDED)
  set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
  find_package (SeqAn3 REQUIRED)
endif()

set(CMAKE_BENCHMARK_ENABLE_GTEST_TESTS OFF)

macro(dependency name)
  # Download and unpack googletest at configure time
  configure_file(CMakeLists.txt.${name} ${name}-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name}-download )
  if(result)
    message(FATAL_ERROR "CMake step for ${name} failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name}-download )
  if(result)
    message(FATAL_ERROR "Build step for ${name} failed: ${result}")
  endif()

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/${name}-src
                   ${CMAKE_CURRENT_BINARY_DIR}/${name}-build
                   EXCLUDE_FROM_ALL)
endmacro(dependency)

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

dependency(googletest)
dependency(benchmark)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

macro(test_case unit_test)
  get_filename_component(target "${unit_test}" NAME_WE)

  add_executable(${target} ${unit_test})
  target_link_libraries(${target} gtest_main)
  target_link_libraries(${target} seqan3::seqan3)
  target_link_libraries(${target} Threads::Threads)

  add_test(NAME ${target} COMMAND ${target})

  unset (unit_test)
  unset (target)
endmacro(test_case)

test_case(probabilistic_suffix_tree_tests.cpp)

add_subdirectory(search)
add_subdirectory(search/lazy_suffix_tree)

add_executable(benchmarks benchmarks.cpp)
target_link_libraries(benchmarks benchmark)
target_link_libraries(benchmarks seqan3::seqan3)
target_link_libraries(benchmarks Threads::Threads)
