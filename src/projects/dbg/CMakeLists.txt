project(debruijn)
set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
add_executable(jumboDB dbg.cpp)
target_link_libraries(jumboDB common sequence m)

