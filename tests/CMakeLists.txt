cmake_minimum_required(VERSION 3.5)


project(uvwmaterial_test)

set(CMAKE_CXX_STANDARD 17)

enable_testing()
include(CTest)

#google test
include(FetchContent)
FetchContent_Declare(
  googletest
#  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
#


include(GoogleTest)
macro(add_uvwmaterial_test TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} uvwmaterial gtest_main)
    add_test(NAME ${TARGET_NAME} COMMAND ${TARGET_NAME})
    #gtest_discover_tests(${TARGET_NAME})
endmacro()

add_uvwmaterial_test(uvwmaterial_test main.cpp)
