/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef THERMAL_MATERIAL_BASE_MEAT_H
#define THERMAL_MATERIAL_BASE_MEAT_H

template<typename _T, std::size_t _Dim, typename _Container>
constexpr thermal_material_base<_T, _Dim, _Container>::thermal_material_base():
    material_base<_T, _Dim, _Container>(),
    _q(),
    _k()
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr thermal_material_base<_T, _Dim, _Container>::thermal_material_base(_Container const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _q(),
    _k()
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename _Parameter>
constexpr thermal_material_base<_T, _Dim, _Container>::thermal_material_base(std::initializer_list<_Parameter> const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _q(),
    _k()
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename ..._Parameter>
constexpr thermal_material_base<_T, _Dim, _Container>::thermal_material_base(_Parameter&& ...__parameter):
    material_base<_T, _Dim, _Container>(__parameter...),
    _q(),
    _k()
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermal_material_base<_T, _Dim, _Container>::heat_flux_tensor()const{
    return _q;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermal_material_base<_T, _Dim, _Container>::heat_flux_tensor(){
    return _q;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermal_material_base<_T, _Dim, _Container>::conductivity_tensor()const{
    return _k;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermal_material_base<_T, _Dim, _Container>::conductivity_tensor(){
    return _k;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermal_material_base<_T, _Dim, _Container>::thermal_gradient()const{
    return _dT;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermal_material_base<_T, _Dim, _Container>::thermal_gradient(){
    return _dT;
}

#endif // THERMAL_MATERIAL_BASE_MEAT_H
