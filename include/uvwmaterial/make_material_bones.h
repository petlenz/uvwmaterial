/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef MAKE_MATERIAL_BONES_H
#define MAKE_MATERIAL_BONES_H

#include "yield_function_base_bones.h"

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_solid_material(material_base<T, Dim, Container> * material){
    return static_cast<solid_material_base<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_solid_material(const material_base<T, Dim, Container> * material){
    return static_cast<const solid_material_base<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_conductivity_material(material_base<T, Dim, Container> * material){
    return static_cast<conductivity_material_base<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_conductivity_material(const material_base<T, Dim, Container> * material){
    return static_cast<const conductivity_material_base<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_small_strain_material(material_base<T, Dim, Container> * material){
    return static_cast<small_strain_material_base<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_small_strain_material(solid_material_base<T, Dim, Container> * material){
    return static_cast<small_strain_material_base<T, Dim, Container>*>(material);
}

//template<typename T, std::size_t Dim, typename Container>
//constexpr inline auto* make_plastic(material_base<T, Dim, Container> * material){
//    return dynamic_cast<plastic_material_base<T, Dim>*>(material);
//}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_history(material_base<T, Dim, Container> * material){
    return dynamic_cast<history_material_base<T>*>(material);
}

//template<typename T, std::size_t Dim, typename Container>
//constexpr inline auto* make_temperature(material_base<T, Dim, Container> * material){
//    return dynamic_cast<temperature_material_base<T>*>(material);
//}

//template<typename T, std::size_t Dim, typename Container>
//constexpr inline auto* make_composite_history(material_base<T, Dim, Container> * material){
//    return dynamic_cast<composite_history_material_base<T, Dim, Container>*>(material);
//}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_yield_function_isotropic_damage_strain_based_base(yield_function_base<T, Dim, Container> * material){
    return dynamic_cast<yield_function_isotropic_damage_strain_based<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_small_strain_material_isotropic_damage(material_base<T, Dim, Container> * material){
    return dynamic_cast<small_strain_isotropic_damage<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_small_strain_material_isotropic_nonlocal_damage(material_base<T, Dim, Container> * material){
    return dynamic_cast<small_strain_isotropic_nonlocal_damage<T, Dim, Container>*>(material);
}

//template<typename T, std::size_t Dim, typename Container>
//constexpr inline auto* make_small_strain_material_eigenstrain_isotropic_nonlocal_damage(material_base<T, Dim, Container> * material){
//    return dynamic_cast<small_strain_eigenstrain_isotropic_nonlocal_damage<T, Dim, Container>*>(material);
//}

//template<typename T, std::size_t Dim, typename Container>
//constexpr inline auto* make_small_strain_material_plasticity_single_yield_function(material_base<T, Dim, Container> * material){
//    return dynamic_cast<small_strain_plasticity_single_yield_function<T, Dim, Container>*>(material);
//}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_yield_function_stress_base(yield_function_base<T, Dim, Container> * material){
    return dynamic_cast<yield_function_stress_base<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto* make_composite_material(material_base<T, Dim, Container> * material){
    return dynamic_cast<composite_material_base<T, Dim, Container>*>(material);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_linear_elasticity(Parameter && ...parameter){
    return std::make_unique<linear_elasticity<T, Dim, Container>>(parameter...);
}

//template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
//static constexpr inline auto make_linear_thermo_elasticity(Parameter && ...parameter){
//    return std::make_unique<linear_thermo_elasticity<T, Dim, Container>>(parameter...);
//}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_small_strain_material_isotropic_damage(Parameter && ...parameter){
    return std::make_unique<small_strain_isotropic_damage<T, Dim, Container>>(parameter...);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_small_strain_material_isotropic_nonlocal_damage(Parameter && ...parameter){
    return std::make_unique<small_strain_isotropic_nonlocal_damage<T, Dim, Container>>(parameter...);
}


//template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
//static constexpr inline auto make_small_strain_material_eigenstrain_isotropic_nonlocal_damage(Parameter && ...parameter){
//    return std::make_unique<small_strain_eigenstrain_isotropic_nonlocal_damage<T, Dim, Container>>(parameter...);
//}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_small_strain_material_plasticity_single_yield_function(Parameter && ...parameter){
    return std::make_unique<small_strain_plasticity_single_yield_function<T, Dim, Container>>(parameter...);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_isotropic_damage_strain_based_yield_function(Parameter && ...parameter){
    return std::make_unique<yield_function_isotropic_damage_strain_based<T, Dim, Container>>(parameter...);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_von_mises_strain_state_function(Parameter && ...parameter){
    return std::make_unique<von_mises_strain_state_function<T, Dim, Container>>(parameter...);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_von_mises_state_function(Parameter && ...parameter){
    return std::make_unique<von_mises_state_function<T, Dim, Container>>(parameter...);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_vector_strain_state_function(Parameter && ...parameter){
    return std::make_unique<vector_strain_state_function<T, Dim, Container>>(parameter...);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_strain_energy_strain_state_function(Parameter && ...parameter){
    return std::make_unique<strain_energy_strain_state_function<T, Dim, Container>>(parameter...);
}

template<typename T, typename ...Parameter>
static constexpr inline auto make_propagation_function_linear(Parameter && ...parameter){
    return std::make_unique<linear<T>>(parameter...);
}

template<typename T, typename ...Parameter>
static constexpr inline auto make_propagation_function_power_law(Parameter && ...parameter){
    return std::make_unique<power_law<T>>(parameter...);
}

template<typename T, typename ...Parameter>
static constexpr inline auto make_propagation_function_exponential_law(Parameter && ...parameter){
    return std::make_unique<exponential_law<T>>(parameter...);
}

template<typename T, typename ...Parameter>
static constexpr inline auto make_propagation_function_exponential1_law(Parameter && ...parameter){
    return std::make_unique<damage_exponential1<T>>(parameter...);
}

template<typename T, typename ...Parameter>
static constexpr inline auto make_propagation_function_exponential2_law(Parameter && ...parameter){
    return std::make_unique<damage_exponential2<T>>(parameter...);
}

template<typename T, typename ...Parameter>
static constexpr inline auto make_propagation_function_exponential3_law(Parameter && ...parameter){
    return std::make_unique<damage_exponential3<T>>(parameter...);
}

template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
static constexpr inline auto make_j2_yield_function(Parameter && ...parameter){
    return std::make_unique<j2_yield_function<T, Dim, Container>>(parameter...);
}

//inclusion
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline std::unique_ptr<inclusion_base<T, Dim, Container>> make_inclusion(Parameter && ...parameter){
    return std::make_unique<inclusion<T, Dim, Container>>(parameter...);
}

// rule_of_mixture
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_rule_of_mixture(Parameter && ...parameter){
    return std::make_unique<small_strain_rule_of_mixture<T, Dim, Container>>(parameter...);
}

// rule_of_mixture history
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_rule_of_mixture_history(Parameter && ...parameter){
    return std::make_unique<small_strain_rule_of_mixture_history<T, Dim, Container>>(parameter...);
}

//rule_of_mixture kernal voigt
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_rule_of_mixture_kernal_voigt(Parameter && ...parameter){
    return std::make_unique<rule_of_mixture_voigt<T, Dim, Container>>(parameter...);
}

//rule_of_mixture kernal reuss
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_rule_of_mixture_kernal_reuss(Parameter && ...parameter){
    return std::make_unique<rule_of_mixture_reuss<T, Dim, Container>>(parameter...);
}

//rule_of_mixture kernal  vrh
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_rule_of_mixture_kernal_vrh(Parameter && ...parameter){
    return std::make_unique<rule_of_mixture_vrh<T, Dim, Container>>(parameter...);
}

//small strain mean field composite
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_mean_field_composite(Parameter && ...parameter){
    return std::make_unique<small_strain_mean_field_composite<T, Dim, Container>>(parameter...);
}

//small strain mean field composite history
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_mean_field_composite_history(Parameter && ...parameter){
    return std::make_unique<small_strain_mean_field_composite_history<T, Dim, Container>>(parameter...);
}

//small strain mean field composite nonlocal damage
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_mean_field_composite_damage_nonlocal(Parameter && ...parameter){
    return std::make_unique<small_strain_mean_field_composite_damage_nonlocal<T, Dim, Container>>(parameter...);
}

//mean field kernal mori tanaka
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_mean_field_mori_tanka_kernal(Parameter && ...parameter){
    return std::make_unique<mean_field_composite_mori_tanaka_kernal<T, Dim, Container>>(parameter...);
}

//mean field kernal dilute
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_mean_field_dilute_kernal(Parameter && ...parameter){
    return std::make_unique<mean_field_composite_dilute_kernal<T, Dim, Container>>(parameter...);
}

//mean field kernal scs
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_mean_field_scs_kernal(Parameter && ...parameter){
    return std::make_unique<mean_field_composite_scs_kernal<T, Dim, Container>>(parameter...);
}

//small strain short fibre composite
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_short_fibre(Parameter && ...parameter){
    return std::make_unique<small_strain_short_fibre_composite<T, Dim, Container>>(parameter...);
}

//small strain short fibre history composite
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_short_fibre_history(Parameter && ...parameter){
    return std::make_unique<small_strain_short_fibre_history_composite<T, Dim, Container>>(parameter...);
}

//small strain short fibre nonlocal damage composite
template<typename T, std::size_t Dim, typename Container, typename ...Parameter>
constexpr inline auto make_small_strain_short_fibre_damage_nonlocal(Parameter && ...parameter){
    return std::make_unique<small_strain_short_fibre_damage_nonlocal_composite<T, Dim, Container>>(parameter...);
}

//eshelby tensor
template<typename T, std::size_t Dim, typename ...Parameter>
static constexpr inline auto make_eshelby_solid_sphere(Parameter && ...parameter){
    return std::make_unique<eshelby_tensor_solid_sphere<T, Dim>>(parameter...);
}

template<typename T, std::size_t Dim, typename ...Parameter>
static constexpr inline auto make_eshelby_long_fibre_x1(Parameter && ...parameter){
    return std::make_unique<eshelby_tensor_solid_cylinder_x1<T, Dim>>(parameter...);
}

template<typename T, std::size_t Dim, typename ...Parameter>
static constexpr inline auto make_eshelby_long_fibre_x2(Parameter && ...parameter){
    return std::make_unique<eshelby_tensor_solid_cylinder_x2<T, Dim>>(parameter...);
}

template<typename T, std::size_t Dim, typename ...Parameter>
static constexpr inline auto make_eshelby_long_fibre_x3(Parameter && ...parameter){
    return std::make_unique<eshelby_tensor_solid_cylinder_x3<T, Dim>>(parameter...);
}

template<typename T, std::size_t Dim, typename ...Parameter>
static constexpr inline auto make_eshelby_long_fibre(Parameter && ...parameter){
    return std::make_unique<eshelby_tensor_solid_cylinder<T, Dim>>(parameter...);
}

template<typename T, std::size_t Dim, typename ...Parameter>
static constexpr inline auto make_eshelby_tensor_solid_cylindrical_ellipsoid(Parameter && ...parameter){
    return std::make_unique<eshelby_tensor_solid_cylindrical_ellipsoid<T, Dim>>(parameter...);
}

#endif // MAKE_MATERIAL_BONES_H
