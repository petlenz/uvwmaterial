/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_CYLINDER_MEAT_H
#define ESHELBY_TENSOR_SOLID_CYLINDER_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder<_T, _Dim>::eshelby_tensor_solid_cylinder():
    _n()
{}

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder<_T, _Dim>::~eshelby_tensor_solid_cylinder(){}

template <typename _T, std::size_t _Dim>
inline void eshelby_tensor_solid_cylinder<_T, _Dim>::init(){
    if(!this->_is_init){
        if constexpr (_Dim == 2){
            tmech::tensor<value_type, 3, 4> Stemp;
            setup_entries(Stemp, this->_parameter[0], _n);
            this->_S = convert_3D_to_2D(Stemp);
        }else{
            setup_entries(this->_S, this->_parameter[0], _n);
        }
        this->_is_init = true;
    }
}

template <typename _T, std::size_t _Dim>
constexpr inline auto& eshelby_tensor_solid_cylinder<_T, _Dim>::direction(){
    return _n;
}

template <typename _T, std::size_t _Dim>
constexpr inline auto const& eshelby_tensor_solid_cylinder<_T, _Dim>::direction()const{
    return _n;
}

template <typename _T, std::size_t _Dim>
template<typename _Tensor, typename _N>
constexpr inline auto eshelby_tensor_solid_cylinder<_T, _Dim>::setup_entries(_Tensor& _S, _T const __nue, _N const& __n){
    const tmech::eye<_T, 3, 2> I;
    const auto p{tmech::otimes(__n,__n)};
    const auto II1{0.5*(tmech::otimesu(I,I) + tmech::otimesl(I,I))};
    const auto II2{tmech::otimes(I,I)};
    const auto II3{tmech::otimes(p,I)};
    const auto II4{tmech::otimes(I,p)};
    const auto II5{(tmech::otimesu(p,I) + tmech::otimesl(p,I) + tmech::otimesu(I,p) + tmech::otimesl(I,p))*0.5};
    const auto II6{tmech::otimes(p, p)};

    const _T c1{(4*__nue-3)/(4*(__nue-1))};
    const _T c2{(4*__nue-1)/(8*(1-__nue))};
    const _T c3{(4*__nue-1)/(8*(__nue-1))};
    const _T c4{1/(8*(1-__nue))};
    const _T c5{(2*__nue-1)/(4*(1-__nue))};
    const _T c6{3/(8*(__nue-1))};
    _S = c1*II1 + c2*II2 + c3*II3 + c4*II4 + c5*II5 + c6*II6;
}


#endif // ESHELBY_TENSOR_SOLID_CYLINDER_MEAT_H
