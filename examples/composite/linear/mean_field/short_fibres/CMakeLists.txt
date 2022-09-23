cmake_minimum_required(VERSION 3.5)

project(short_fibres LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)

add_executable(short_fibres main.cpp)
target_link_libraries(short_fibres PUBLIC matplot)
