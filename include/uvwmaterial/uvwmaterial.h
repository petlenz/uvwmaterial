#ifndef UVWMATERIAL_H
#define UVWMATERIAL_H

#include <memory>
#include <cstdlib>
#include <vector>
#include <map>

#include <tmech/tmech.h>

#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>

namespace uvwmat {


template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp){
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
            // unless the result is subnormal
            || std::fabs(x-y) < std::numeric_limits<T>::min();
}

namespace detail {
#include "uvwmat_utility.h"
}

//lebedev integration points
#include "sphere_lebedev_rule_bones.h"
#include "sphere_lebedev_rule_meat.h"

//functions bones
#include "uvwmat_functions_bones.h"

//base material
#include "material_base_bones.h"
#include "material_base_meat.h"

#include "solid_material_base_bones.h"
#include "solid_material_base_meat.h"
#include "conductivity_material_base_bones.h"
#include "conductivity_material_base_meat.h"
#include "thermal_material_base_bones.h"
#include "thermal_material_base_meat.h"
#include "fluid_material_base_bones.h"
#include "fluid_material_base_meat.h"
#include "electro_magnetic_material_base_bones.h"
#include "electro_magnetic_material_base_meat.h"


//curing functions
#include "polymer_curing_functions_bones.h"
#include "polymer_curing_functions_meat.h"

//incremental material
#include "incremental_solid_material_bones.h"
#include "incremental_solid_material_meat.h"

//time dependent material
#include "time_dependent_material_base_bones.h"
#include "time_dependent_material_base_meat.h"

//temperature dependent material
#include "temperature_dependent_material_base_bones.h"
#include "temperature_dependent_material_base_meat.h"

//degree of cure dependent material
#include "degree_of_cure_dependent_material_bones.h"
#include "degree_of_cure_dependent_material_meat.h"

//history material base
#include "history_material_base_bones.h"
#include "history_material_base_meat.h"

//plastic material base
#include "plastic_material_base_bones.h"
#include "plastic_material_base_meat.h"

//nonlocal material
#include "nonlocal_material_base_bones.h"
#include "nonlocal_material_base_meat.h"

//eigenstrain material
#include "eigenstrain_material_base_bones.h"
#include "eigenstrain_material_base_meat.h"


#include "thermo_mechanical_material_base_bones.h"
#include "thermo_mechanical_material_base_meat.h"

//finite strain material
#include "finite_strain_solid_material_base_bones.h"
#include "finite_strain_solid_material_base_meat.h"

//eigenstrain material
#include "eigenstrain_material_base_bones.h"
#include "eigenstrain_material_base_meat.h"
#include "eigenstrain_temperature_material_base_bones.h"
#include "eigenstrain_temperature_material_base_meat.h"
#include "eigenstrain_curing_material_base_bones.h"
#include "eigenstrain_curing_material_base_meat.h"

//state functions
#include "state_function_base_bones.h"
#include "state_function_base_meat.h"
#include "state_function_strain_based_bones.h"
#include "state_function_strain_based_meat.h"
#include "state_function_stress_based_bones.h"
#include "state_function_stress_based_meat.h"

//propagation laws
#include "propagation_law_base_bones.h"
#include "propagation_law_base_meat.h"
#include "propagation_strain_law_bones.h"
#include "propagation_strain_law_meat.h"
#include "propagation_law_viscous_regularization_bones.h"
#include "propagation_law_viscous_regularization_meat.h"

//yield functions
#include "yield_function_base_bones.h"
#include "yield_function_base_meat.h"

//tensor isotropization
#include "tensor_isotropization_base_bones.h"
#include "tensor_isotropization_base_meat.h"
#include "tensor_isotropization_general_bones.h"
#include "tensor_isotropization_general_meat.h"


//eshelby tensor
#include "eshelby_tensor_functions.h"
#include "eshelby_tensor_base_bones.h"
#include "eshelby_tensor_base_meat.h"
#include "eshelby_tensor_solid_base_bones.h"
#include "eshelby_tensor_solid_base_meat.h"
#include "eshelby_tensor_conductivity_base_bones.h"
#include "eshelby_tensor_conductivity_base_meat.h"
#include "eshelby_tensor_nonlinear_base_bones.h"
#include "eshelby_tensor_nonlinear_base_meat.h"
#include "eshelby_tensor_solid_cylindrical_ellipsoid_bones.h"
#include "eshelby_tensor_solid_cylindrical_ellipsoid_meat.h"
#include "eshelby_tensor_solid_cylinder_bones.h"
#include "eshelby_tensor_solid_cylinder_meat.h"
#include "eshelby_tensor_solid_cylinder_x1_bones.h"
#include "eshelby_tensor_solid_cylinder_x1_meat.h"
#include "eshelby_tensor_solid_cylinder_x2_bones.h"
#include "eshelby_tensor_solid_cylinder_x2_meat.h"
#include "eshelby_tensor_solid_cylinder_x3_bones.h"
#include "eshelby_tensor_solid_cylinder_x3_meat.h"
#include "eshelby_tensor_solid_sphere_bones.h"
#include "eshelby_tensor_solid_sphere_meat.h"
#include "eshelby_tensor_solid_general_geometry_bones.h"
#include "eshelby_tensor_solid_general_geometry_meat.h"
#include "eshelby_tensor_solid_sphere_nonlinear_small_strain_bones.h"
#include "eshelby_tensor_solid_sphere_nonlinear_small_strain_meat.h"
#include "eshelby_tensor_conductivity_sphere_bones.h"
#include "eshelby_tensor_conductivity_sphere_meat.h"
#include "eshelby_tensor_finite_strain_base_bones.h"
#include "eshelby_tensor_finite_strain_base_meat.h"
#include "eshelby_tensor_solid_general_geometry_finite_strain_bones.h"
#include "eshelby_tensor_solid_general_geometry_finite_strain_meat.h"

//inclusion base
#include "composite_material_inclusion_base_bones.h"
#include "composite_material_inclusion_base_meat.h"
#include "inclusion_bones.h"
#include "inclusion_meat.h"

//composite material base
#include "composite_material_base_bones.h"
#include "composite_material_base_meat.h"
#include "composite_material_solid_base_bones.h"
#include "composite_material_solid_base_meat.h"
#include "composite_material_conductivity_base_bones.h"
#include "composite_material_conductivity_base_meat.h"
#include "composite_material_history_base_bones.h"
#include "composite_material_history_base_meat.h"
#include "composite_material_nonlinear_base_bones.h"
#include "composite_material_nonlinear_base_meat.h"
#include "composite_material_eigenstrain_base_bones.h"
#include "composite_material_eigenstrain_base_meat.h"
#include "nonlocal_composite_base_bones.h"
#include "nonlocal_composite_base_meat.h"

//short fibre material
#include "composite_material_short_fibre_material_base_bones.h"
#include "composite_material_short_fibre_material_base_meat.h"
#include "short_fibre_history_base_bones.h"
#include "short_fibre_history_base_meat.h"
#include "short_fibre_nonlocal_base_bones.h"
#include "short_fibre_nonlocal_base_meat.h"

// small strain material base
#include "small_strain_material_base_bones.h"
#include "small_strain_material_base_meat.h"

//mean field kernal
#include "mean_field_composite_kernal_base_bones.h"
#include "mean_field_composite_kernal_base_meat.h"
#include "mean_field_composite_solid_kernal_base_bones.h"
#include "mean_field_composite_solid_kernal_base_meat.h"
#include "mean_field_composite_solid_dilute_kernal_bones.h"
#include "mean_field_composite_solid_dilute_kernal_meat.h"
#include "mean_field_composite_solid_mori_tanaka_kernal_bones.h"
#include "mean_field_composite_solid_mori_tanaka_kernal_meat.h"
#include "mean_field_composite_solid_scs_kernal_bones.h"
#include "mean_field_composite_solid_scs_kernal_meat.h"
#include "mean_field_composite_base_bones.h"
#include "mean_field_composite_base_meat.h"
#include "mean_field_composite_conductivity_kernal_base_bones.h"
#include "mean_field_composite_conductivity_kernal_base_meat.h"
#include "mean_field_composite_conductivity_base_bones.h"
#include "mean_field_composite_conductivity_base_meat.h"
#include "mean_field_composite_solid_base_bones.h"
#include "mean_field_composite_solid_base_meat.h"
#include "mean_field_composite_conductivity_material_bones.h"
#include "mean_field_composite_conductivity_material_meat.h"

#include "mean_field_composite_solid_finite_strain_dilute_kernal_bones.h"
#include "mean_field_composite_solid_finite_strain_dilute_kernal_meat.h"

//rule of mixture kernal
#include "composite_material_rule_of_mixture_kernal_base_bones.h"
#include "composite_material_rule_of_mixture_kernal_base_meat.h"
#include "composite_material_rule_of_mixture_base_bones.h"
#include "composite_material_rule_of_mixture_base_meat.h"

//small strain materials
#include "linear_elasticity_bones.h"
#include "linear_elasticity_meat.h"
#include "small_strain_mean_field_material_bones.h"
#include "small_strain_mean_field_material_meat.h"
#include "small_strain_rule_of_mixture_material_bones.h"
#include "small_strain_rule_of_mixture_material_meat.h"
#include "small_strain_short_fibre_composite_material_bones.h"
#include "small_strain_short_fibre_composite_material_meat.h"
#include "small_strain_short_fibre_history_composite_bones.h"
#include "small_strain_short_fibre_history_composite_meat.h"
#include "small_strain_short_fibre_nonlocal_composite_bones.h"
#include "small_strain_short_fibre_nonlocal_composite_meat.h"
#include "small_strain_isotropic_damage_material_bones.h"
#include "small_strain_isotropic_damage_material_meat.h"
#include "small_strain_rule_of_mixture_history_base_bones.h"
#include "small_strain_rule_of_mixture_history_base_meat.h"
#include "small_strain_mean_field_history_bones.h"
#include "small_strain_mean_field_history_meat.h"
#include "small_strain_isotropic_nonlocal_damage_bones.h"
#include "small_strain_isotropic_nonlocal_damage_meat.h"
#include "small_strain_mean_field_nonlocal_bones.h"
#include "small_strain_mean_field_nonlocal_meat.h"
#include "small_strain_plasticity_single_yield_function_bones.h"
#include "small_strain_plasticity_single_yield_function_meat.h"

//conductivity
#include "polymer_matrix_conductivity_material_bones.h"


//incremental
#include "incremental_linear_elasticity_bones.h"
#include "incremental_linear_elasticity_meat.h"


//thermal material
#include "linear_isotropic_thermal_conductivity_bones.h"
#include "linear_isotropic_thermal_conductivity_meat.h"

//finite strain materials
#include "finite_strain_mean_field_composite_bones.h"
#include "finite_strain_mean_field_composite_meat.h"

#include "saint_venant_kirchhoff_bones.h"
#include "saint_venant_kirchhoff_meat.h"
#include "incremental_saint_venant_kirchhoff_bones.h"
#include "incremental_saint_venant_kirchhoff_meat.h"
#include "hencky_strains_linear_elastic_bones.h"
#include "hencky_strains_linear_elastic_bones.h"


//plastic material
//#include "plastic_material_base_bones.h"
//#include "plastic_material_base_meat.h"
//#include "plastic_material_anisotropic_yielding_base_bones.h"
//#include "plastic_material_anisotropic_yielding_base_meat.h"
//#include "plastic_material_isotropic_yielding_base_bones.h"
//#include "plastic_material_isotropic_yielding_base_meat.h"

//make unique ptr
#include "make_material_bones.h"

//

#include "thermo_mechanical_elastic_bones.h"
#include "thermo_mechanical_elastic_meat.h"

#include "gamm2022_bones.h"
#include "gamm2022_meat.h"



//general material
#include "general_material_bones.h"
#include "general_material_meat.h"

//functions meat
#include "uvwmat_functions_meat.h"

//element deletion
#include "element_deletion_bones.h"
#include "element_deletion_meat.h"



}

#endif // UVWMATERIAL_H
