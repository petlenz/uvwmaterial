cmake_minimum_required(VERSION 3.5)

project(mori_tanaka LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Matplot++ REQUIRED)

add_executable(mori_tanaka main.cpp)
target_link_libraries(mori_tanaka PUBLIC matplot)
