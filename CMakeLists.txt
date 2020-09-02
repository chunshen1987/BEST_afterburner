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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -g -O3 -std=c++11 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

add_executable(sampler_and_smash src/SamplerAndSmash.cc)
include_directories(src/include)

message("SMASH_DIR = $ENV{SMASH_DIR}")
find_package(SMASH)
if(${SMASH_FOUND})
    include_directories(SYSTEM ${SMASH_INCLUDE_DIR})
    include_directories(SYSTEM $ENV{SMASH_DIR}/build/src/include)
    target_link_libraries(sampler_and_smash ${SMASH_LIBRARIES})
endif(${SMASH_FOUND})

option(iSS "Build tests for iSS" OFF)

set(SAMPLER_INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}/external_codes/microcanonical_cooper_frye/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external_codes/best_sampler/software/include)

if (iSS)
    set(SAMPLER_INCLUDES
        ${SAMPLER_INCLUDES}
        ${CMAKE_CURRENT_SOURCE_DIR}/external_codes/iSS
        ${CMAKE_CURRENT_SOURCE_DIR}/external_codes/iSS/src)
    set_target_properties (sampler_and_smash PROPERTIES COMPILE_FLAGS "-DiSSFlag")
endif (iSS)

include_directories(${SAMPLER_INCLUDES})
message("SAMPLER_INCLUDES:  ${SAMPLER_INCLUDES}")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external_codes/microcanonical_cooper_frye")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external_codes/best_sampler")

target_link_libraries(sampler_and_smash microcanonical_sampling msu_sampler_shared)

if (iSS)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external_codes/iSS")
    target_link_libraries(sampler_and_smash iSS)
endif (iSS)

option(USE_ROOT "Turn this off to disable ROOT output support in SMASH." OFF)
if (USE_ROOT)
    find_package(ROOT 5.34)
    if (ROOT_FOUND)
        include_directories(SYSTEM "${ROOT_INCLUDE_DIR}")
        add_definitions(-DSMASH_USE_ROOT)
        target_link_libraries(sampler_and_smash ${ROOT_LIBRARIES})
    endif ()
endif (USE_ROOT)
