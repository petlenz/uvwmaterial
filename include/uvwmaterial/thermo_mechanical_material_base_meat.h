/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef THERMO_MECHANICAL_MATERIAL_BASE_MEAT_H
#define THERMO_MECHANICAL_MATERIAL_BASE_MEAT_H


template<typename _T, std::size_t _Dim, typename _Container>
constexpr thermo_mechanical_material_base<_T, _Dim, _Container>::thermo_mechanical_material_base():
    material_base<_T, _Dim, _Container>(),
    _stress(),
    _C(),
    _q(),
    _k(),
    _cd(0),
    _rho(0),
    _temperature(0),
    _delta_temperature(0),
    _latent_heat(0),
    _dlatent_heat(0)
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr thermo_mechanical_material_base<_T, _Dim, _Container>::thermo_mechanical_material_base(_Container const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _stress(),
    _C(),
    _q(),
    _k(),
    _cd(0),
    _rho(0),
    _temperature(0),
    _delta_temperature(0),
    _latent_heat(0),
    _dlatent_heat(0)
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename _Parameter>
constexpr thermo_mechanical_material_base<_T, _Dim, _Container>::thermo_mechanical_material_base(std::initializer_list<_Parameter> const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _stress(),
    _C(),
    _q(),
    _k(),
    _cd(0),
    _rho(0),
    _temperature(0),
    _delta_temperature(0),
    _latent_heat(0),
    _dlatent_heat(0)
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename ..._Parameter>
constexpr thermo_mechanical_material_base<_T, _Dim, _Container>::thermo_mechanical_material_base(_Parameter&& ...__parameter):
    material_base<_T, _Dim, _Container>(__parameter...),
    _stress(),
    _C(),
    _q(),
    _k(),
    _cd(0),
    _rho(0),
    _temperature(0),
    _delta_temperature(0),
    _latent_heat(0)
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermo_mechanical_material_base<_T, _Dim, _Container>::heat_flux_tensor()const{
    return _q;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermo_mechanical_material_base<_T, _Dim, _Container>::heat_flux_tensor(){
    return _q;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermo_mechanical_material_base<_T, _Dim, _Container>::conductivity_tensor()const{
    return _k;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermo_mechanical_material_base<_T, _Dim, _Container>::conductivity_tensor(){
    return _k;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermo_mechanical_material_base<_T, _Dim, _Container>::thermal_gradient()const{
    return _dT;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermo_mechanical_material_base<_T, _Dim, _Container>::thermal_gradient(){
    return _dT;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermo_mechanical_material_base<_T, _Dim, _Container>::stress_tensor()const{
    return _stress;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermo_mechanical_material_base<_T, _Dim, _Container>::stress_tensor(){
    return _stress;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermo_mechanical_material_base<_T, _Dim, _Container>::tangent_tensor()const{
    return _C;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermo_mechanical_material_base<_T, _Dim, _Container>::tangent_tensor(){
    return _C;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermo_mechanical_material_base<_T, _Dim, _Container>::temperature(){
    return _temperature;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermo_mechanical_material_base<_T, _Dim, _Container>::temperature()const{
    return _temperature;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& thermo_mechanical_material_base<_T, _Dim, _Container>::delta_temperature(){
    return _delta_temperature;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& thermo_mechanical_material_base<_T, _Dim, _Container>::delta_temperature()const{
    return _delta_temperature;
}
#endif // THERMO_MECHANICAL_MATERIAL_BASE_MEAT_H
