cmake_minimum_required(VERSION 3.5)

project(vrh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)

add_executable(vrh main.cpp)
target_link_libraries(vrh PUBLIC matplot)
