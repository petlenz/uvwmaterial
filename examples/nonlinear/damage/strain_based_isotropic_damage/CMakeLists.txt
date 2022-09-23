cmake_minimum_required(VERSION 3.5)

project(strain_based_damage LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)

add_executable(strain_based_damage main.cpp)
target_link_libraries(strain_based_damage PUBLIC matplot)
