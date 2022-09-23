cmake_minimum_required(VERSION 3.11)


# Project
# =======
project(uvwmaterial CXX)

set(UVWMATERIAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(UVWMATERIAL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/uvwmaterial)


# Dependencies
# ============
find_package(tmech 1.0.0 REQUIRED)
find_package(Boost 1.71.0 REQUIRED)
find_package(nlohmann_json 3.1.1 REQUIRED)

# Versionning
# ===========

file(STRINGS "${UVWMATERIAL_INCLUDE_DIR}/uvwmaterial_config.h" uvwmaterial_version_defines
     REGEX "#define UVWMATERIAL_VERSION_(MAJOR|MINOR|PATCH)")
foreach(ver ${uvwmaterial_version_defines})
    if(ver MATCHES "#define UVWMATERIAL_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(UVWMATERIAL_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()
set(${PROJECT_NAME}_VERSION
    ${UVWMATERIAL_VERSION_MAJOR}.${UVWMATERIAL_VERSION_MINOR}.${UVWMATERIAL_VERSION_PATCH})
message(STATUS "Building uvwmaterial v${${PROJECT_NAME}_VERSION}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#tensor calculus
find_package(tmech REQUIRED)


set(UVWMATERIAL_HEADER
    ${UVWMATERIAL_INCLUDE_DIR}/uvwmat_utility.h
    ${UVWMATERIAL_INCLUDE_DIR}/uvwmat_functions_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/uvwmat_functions_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/general_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/general_material_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_mean_field_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_mean_field_material_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_rule_of_mixture_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_rule_of_mixture_material_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eigenstrain_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eigenstrain_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eigenstrain_temperature_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eigenstrain_temperature_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eigenstrain_curing_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eigenstrain_curing_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/history_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/history_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_yield_function_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_yield_function_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_isotropic_yielding_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_isotropic_yielding_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_anisotropic_yielding_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/plastic_material_anisotropic_yielding_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/sphere_lebedev_rule_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/sphere_lebedev_rule_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_short_fibre_composite_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_short_fibre_composite_material_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/make_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/propagation_law_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/propagation_law_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/propagation_law_viscous_regularization_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/propagation_law_viscous_regularization_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/propagation_strain_law_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/propagation_strain_law_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/yield_function_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/yield_function_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/state_function_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/state_function_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/state_function_strain_based_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/state_function_strain_based_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/state_function_stress_based_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/state_function_stress_based_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_isotropic_damage_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_isotropic_damage_material_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/rule_of_mixture_history_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/rule_of_mixture_history_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_rule_of_mixture_history_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_rule_of_mixture_history_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_mean_field_history_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_mean_field_history_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_mean_field_nonlocal_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_mean_field_nonlocal_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/nonlocal_composite_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/nonlocal_composite_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/nonlocal_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/nonlocal_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_isotropic_nonlocal_damage_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_isotropic_nonlocal_damage_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/time_dependent_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/time_dependent_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/inclusion_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/inclusion_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/element_deletion_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/element_deletion_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/short_fibre_history_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/short_fibre_history_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/short_fibre_nonlocal_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/short_fibre_nonlocal_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_short_fibre_history_composite_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_short_fibre_history_composite_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_short_fibre_nonlocal_composite_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_short_fibre_nonlocal_composite_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/solid_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/solid_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/thermal_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/thermal_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/fluid_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/fluid_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/electro_magnetic_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/electro_magnetic_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/linear_isotropic_thermal_conductivity_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/linear_isotropic_thermal_conductivity_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/incremental_solid_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/incremental_solid_material_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_plasticity_single_yield_function_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/small_strain_plasticity_single_yield_function_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/polymer_curing_functions_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/polymer_curing_functions_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/temperature_dependent_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/temperature_dependent_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/incremental_linear_elasticity_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/incremental_linear_elasticity_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/tensor_isotropization_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/tensor_isotropization_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/tensor_isotropization_general_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/tensor_isotropization_general_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_nonlinear_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_nonlinear_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_x1_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_x1_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_x2_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_x2_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_x3_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_x3_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_functions.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_sphere_nonlinear_small_strain_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_sphere_nonlinear_small_strain_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_finite_strain_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_finite_strain_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_general_geometry_finite_strain_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_general_geometry_finite_strain_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_anisotropic_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_anisotropic_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylindrical_ellipsoid_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylindrical_ellipsoid_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_general_geometry_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_general_geometry_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_conductivity_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_conductivity_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_cylinder_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_sphere_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_solid_sphere_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_conductivity_sphere_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/eshelby_tensor_conductivity_sphere_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_solid_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_solid_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_conductivity_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_conductivity_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_nonlinear_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_nonlinear_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_history_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_history_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_eigenstrain_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_eigenstrain_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_rule_of_mixture_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_rule_of_mixture_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_rule_of_mixture_kernal_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_rule_of_mixture_kernal_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_inclusion_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_inclusion_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_short_fibre_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/composite_material_short_fibre_material_base_meat.h


    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_kernal_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_kernal_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_kernal_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_kernal_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_dilute_kernal_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_dilute_kernal_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_mori_tanaka_kernal_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_mori_tanaka_kernal_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_scs_kernal_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_scs_kernal_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_conductivity_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_conductivity_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_conductivity_kernal_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_conductivity_kernal_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_conductivity_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_conductivity_material_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_finite_strain_dilute_kernal_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/mean_field_composite_solid_finite_strain_dilute_kernal_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/conductivity_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/conductivity_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/polymer_matrix_conductivity_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/polymer_matrix_conductivity_material_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/linear_elasticity_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/linear_elasticity_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/neo_hookean_incompressible_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/neo_hookean_incompressible_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/finite_strain_solid_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/finite_strain_solid_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/finite_strain_mean_field_composite_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/finite_strain_mean_field_composite_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/saint_venant_kirchhoff_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/saint_venant_kirchhoff_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/incremental_saint_venant_kirchhoff_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/incremental_saint_venant_kirchhoff_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/hencky_strains_linear_elastic_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/hencky_strains_linear_elastic_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/thermo_mechanical_material_base_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/thermo_mechanical_material_base_meat.h
    ${UVWMATERIAL_INCLUDE_DIR}/thermo_mechanical_elastic_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/thermo_mechanical_elastic_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/gamm2022_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/gamm2022_meat.h

    ${UVWMATERIAL_INCLUDE_DIR}/degree_of_cure_dependent_material_bones.h
    ${UVWMATERIAL_INCLUDE_DIR}/degree_of_cure_dependent_material_meat.h
   )


# add the library that runs
if(${CMAKE_VERSION} VERSION_LESS "3.19.0")
add_library(${PROJECT_NAME} INTERFACE)
else()
add_library(${PROJECT_NAME} INTERFACE
    ${UVWMATERIAL_INCLUDE_DIR}/uvwmaterial.h
    ${UVWMATERIAL_INCLUDE_DIR}/uvwmaterial_config.h
    ${UVWMATERIAL_HEADER}
    )
endif()
# add alias so the project can be uses with add_subdirectory
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

include(GNUInstallDirs)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(IS_TOPLEVEL_PROJECT TRUE)
else()
  set(IS_TOPLEVEL_PROJECT FALSE)
endif()


OPTION(UVWMATERIAL_INSTALL_LIBRARY
       "Enable installing of uvwmaterial library into default locations"
       ${IS_TOPLEVEL_PROJECT})
OPTION(BUILD_EXAMPLES "uvwmaterial examples" ON)
OPTION(BUILD_TESTS "uvwmaterial test suite" ON)
OPTION(BUILD_BENCHMARK "uvwmaterial benchmark" OFF)
OPTION(DOWNLOAD_GTEST "build gtest from downloaded sources" OFF)
OPTION(DOWNLOAD_GBENCHMARK "download google benchmark and build from source" OFF)



# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUIcd build && ctestnstallDirs)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)


if(UVWMATERIAL_INSTALL_LIBRARY)
  # locations are provided by GNUInstallDirs
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    "${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${${PROJECT_NAME}_VERSION}
    COMPATIBILITY SameMajorVersion)

  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(
    EXPORT ${PROJECT_NAME}_Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
          DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/uvwmaterial
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()


set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}")

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

#if(BUILD_BENCHMARK)
#    add_subdirectory(benchmark)
#endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

include(CPack)
