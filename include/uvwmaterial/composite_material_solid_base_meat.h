/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef COMPOSITE_MATERIAL_SOLID_BASE_MEAT_H
#define COMPOSITE_MATERIAL_SOLID_BASE_MEAT_H


template <typename _T, std::size_t _Dim, typename _Container>
composite_material_solid_base<_T, _Dim, _Container>::composite_material_solid_base(){}

template <typename _T, std::size_t _Dim, typename _Container>
composite_material_solid_base<_T, _Dim, _Container>::~composite_material_solid_base(){}


template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& composite_material_solid_base<_T, _Dim, _Container>::strain_concentration_tensors()const{
    return _strain_concentration_tensors;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& composite_material_solid_base<_T, _Dim, _Container>::strain_concentration_tensors(){
    return _strain_concentration_tensors;
}

#endif // COMPOSITE_MATERIAL_SOLID_BASE_MEAT_H
