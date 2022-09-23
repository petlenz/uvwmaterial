cmake_minimum_required(VERSION 3.5)

project(strain_based_damage_mean_field_composite LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)
find_package(tmech REQUIRED)

add_executable(strain_based_damage_mean_field_composite main.cpp)
target_link_libraries(strain_based_damage_mean_field_composite PUBLIC matplot)
