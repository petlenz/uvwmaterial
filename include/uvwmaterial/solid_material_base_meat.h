/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SOLID_MATERIAL_BASE_MEAT_H
#define SOLID_MATERIAL_BASE_MEAT_H

template<typename _T, std::size_t _Dim, typename _Container>
constexpr solid_material_base<_T, _Dim, _Container>::solid_material_base():
    material_base<_T, _Dim, _Container>(),
    _stress(),
    _C()
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr solid_material_base<_T, _Dim, _Container>::solid_material_base(_Container const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _stress(),
    _C()
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename _Parameter>
constexpr solid_material_base<_T, _Dim, _Container>::solid_material_base(std::initializer_list<_Parameter> const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _stress(),
    _C()
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename ..._Parameter>
constexpr solid_material_base<_T, _Dim, _Container>::solid_material_base(_Parameter&& ...__parameter):
    material_base<_T, _Dim, _Container>(__parameter...),
    _stress(),
    _C()
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& solid_material_base<_T, _Dim, _Container>::stress_tensor()const{
    return _stress;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& solid_material_base<_T, _Dim, _Container>::stress_tensor(){
    return _stress;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& solid_material_base<_T, _Dim, _Container>::tangent_tensor()const{
    return _C;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& solid_material_base<_T, _Dim, _Container>::tangent_tensor(){
    return _C;
}

#endif // SOLID_MATERIAL_BASE_MEAT_H
