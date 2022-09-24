/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef TEMPERATURE_DEPENDENT_MATERIAL_BASE_MEAT_H
#define TEMPERATURE_DEPENDENT_MATERIAL_BASE_MEAT_H

template <typename T>
temperature_dependent_material_base<T>::temperature_dependent_material_base():
    _temperature(0),
    _delta_temperature(0)
{}

template <typename T>
temperature_dependent_material_base<T>::~temperature_dependent_material_base(){}


template <typename T>
constexpr inline auto& temperature_dependent_material_base<T>::temperature(){
    return _temperature;
}

template <typename T>
constexpr inline auto const& temperature_dependent_material_base<T>::temperature()const{
    return _temperature;
}

template <typename T>
constexpr inline auto& temperature_dependent_material_base<T>::delta_temperature(){
    return _delta_temperature;
}

template <typename T>
constexpr inline auto const& temperature_dependent_material_base<T>::delta_temperature()const{
    return _delta_temperature;
}

#endif // TEMPERATURE_DEPENDENT_MATERIAL_BASE_MEAT_H
