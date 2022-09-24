cmake_minimum_required(VERSION 3.5)

project(json_linear_elasticity LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)


add_executable(json_linear_elasticity main.cpp)
target_link_libraries(json_linear_elasticity PUBLIC Matplot++::matplot)

