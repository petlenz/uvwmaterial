cmake_minimum_required(VERSION 3.5)

project(linear_elasticity LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(linear_elasticity main.cpp)

find_package(Matplot++ REQUIRED)

target_link_libraries(linear_elasticity PUBLIC matplot)
