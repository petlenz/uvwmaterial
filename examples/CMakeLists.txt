macro(add_uvwmaterial_example TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    target_link_libraries(${TARGET_NAME} PRIVATE uvwmaterial)
    #maybe_target_pedantic_warnings(${TARGET_NAME})
endmacro()

#add_subdirectory(dir)
