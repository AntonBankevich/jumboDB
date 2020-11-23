cmake_minimum_required(VERSION 3.1)
project(DR)
set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "-std=c++14 -lstdc++fs -ggdb3 ${OpenMP_CXX_FLAGS}" )
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_SHARED_LINKER_FLAGS "-Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -ltbb -fsigned-char")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/tools)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/lib)
find_package(OpenMP)
add_subdirectory(src/tools)
add_subdirectory(src/projects)


