cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (SamplerAndSmash CXX C)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message("source directory:" ${CMAKE_SOURCE_DIR})
  message("binary directory:" ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

set(ENV{SMASH_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/external_codes/smash")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external_codes/smash/cmake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

add_executable(sampler_and_smash src/SamplerAndSmash.cc)
include_directories(src/include)

message($ENV{SMASH_DIR})
find_package(SMASH)
if(${SMASH_FOUND})
  include_directories(SYSTEM ${SMASH_INCLUDE_DIR})
  include_directories(SYSTEM $ENV{SMASH_DIR}/build/src/include)
  target_link_libraries(sampler_and_smash ${SMASH_LIBRARIES})
endif(${SMASH_FOUND})

set(SAMPLER_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/external_codes/microcanonical_cooper_frye/include)
message("SAMPLER_INCLUDES:  ${SAMPLER_INCLUDES}")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external_codes/microcanonical_cooper_frye")
include_directories(${SAMPLER_INCLUDES})

