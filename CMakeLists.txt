cmake_minimum_required(VERSION 3.10) # Establece la versión mínima requerida de CMake
project(RESCOL) # Nombre del proyecto

# Establece el estándar de C++ deseado sin requerirlo explícitamente
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Añade los directorios de los headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/argparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/rapidjson)

# Compila todos los archivos .cpp encontrados en el directorio del proyecto y los enlaza en el ejecutable RESCOL
file(GLOB project_sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_executable(RESCOL ${project_sources})
target_link_libraries(RESCOL stdc++fs pthread)
