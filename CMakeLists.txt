
cmake_minimum_required(VERSION 3.10)
project(RESCOL)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/argparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/rapidjson)

file(GLOB project_sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_executable(RESCOL ${project_sources})

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(RESCOL PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
        -g
    )
endif()

target_link_libraries(RESCOL stdc++fs pthread)