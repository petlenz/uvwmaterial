cmake_minimum_required(VERSION 3.5)

project(reuss LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)

add_executable(reuss main.cpp)
target_link_libraries(reuss PUBLIC matplot)
