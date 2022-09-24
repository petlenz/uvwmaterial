/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef TIME_DEPENDENT_MATERIAL_BASE_MEAT_H
#define TIME_DEPENDENT_MATERIAL_BASE_MEAT_H

template <typename T>
time_dependent_material_base<T>::time_dependent_material_base():
    _time(0),
    _delta_time(0)
{}

template <typename T>
time_dependent_material_base<T>::~time_dependent_material_base(){}


template <typename T>
constexpr inline auto& time_dependent_material_base<T>::time(){
    return _time;
}

template <typename T>
constexpr inline auto const& time_dependent_material_base<T>::time()const{
    return _time;
}

template <typename T>
constexpr inline auto& time_dependent_material_base<T>::delta_time(){
    return _delta_time;
}

template <typename T>
constexpr inline auto const& time_dependent_material_base<T>::delta_time()const{
    return _delta_time;
}

#endif // TIME_DEPENDENT_MATERIAL_BASE_MEAT_H
