/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_CYLINDRICAL_ELLIPSOID_MEAT_H
#define ESHELBY_TENSOR_SOLID_CYLINDRICAL_ELLIPSOID_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylindrical_ellipsoid<_T, _Dim>::eshelby_tensor_solid_cylindrical_ellipsoid(){}

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylindrical_ellipsoid<_T, _Dim>::~eshelby_tensor_solid_cylindrical_ellipsoid(){}

template <typename _T, std::size_t _Dim>
inline void eshelby_tensor_solid_cylindrical_ellipsoid<_T, _Dim>::init(){
    if(!this->_is_init){
        if constexpr (_Dim == 2){
            //std::cout<<_n<<std::endl;
            const tmech::tensor<_T, 3, 4> Stemp{ellipsoid_direction(this->_parameter[0], this->_parameter[1], _n)};
            this->_S = convert_3D_to_2D(Stemp);
        }else{
            //std::cout<<_n<<std::endl;
            this->_S = ellipsoid_direction(this->_parameter[0], this->_parameter[1], _n);
        }
        this->_is_init = true;
    }
}

template <typename _T, std::size_t _Dim>
constexpr inline auto& eshelby_tensor_solid_cylindrical_ellipsoid<_T, _Dim>::direction(){
    return _n;
}

template <typename _T, std::size_t _Dim>
constexpr inline auto const& eshelby_tensor_solid_cylindrical_ellipsoid<_T, _Dim>::direction()const{
    return _n;
}

template <typename _T, std::size_t _Dim>
template<typename _N>
constexpr inline auto eshelby_tensor_solid_cylindrical_ellipsoid<_T, _Dim>::ellipsoid_direction(_T const __nue, _T const __eta, _N const& __n){
    const _T eta2{__eta*__eta};
    //eta = length/diameter
    auto Tg{[__eta, eta2](){
            if(__eta < 1.0){
                return (__eta/(std::pow((1.0-eta2),1.5)))*(std::acos(__eta) - __eta*std::sqrt(1.0-eta2));
            }else{
                return (__eta/(std::pow(eta2-1.0, 1.5)))*(__eta*std::sqrt(eta2 - 1.0) - (std::acosh(__eta)));
            }
        }};

    const _T g{Tg()};

    const tmech::eye<_T, 3, 2> I;
    const auto p{tmech::otimes(__n,__n)};
    const auto II1{0.5*(tmech::otimesu(I,I) + tmech::otimesl(I,I))};
    const auto II2{tmech::otimes(I,I)};
    const auto II3{tmech::otimes(p,I)};
    const auto II4{tmech::otimes(I,p)};
    const auto II5{(tmech::otimesu(p,I) + tmech::otimesl(p,I) + tmech::otimesu(I,p) + tmech::otimesl(I,p))*0.5};
    const auto II6{tmech::otimes(p, p)};

    const _T c{(eta2 - 1)*__nue - eta2 + 1};
    const _T c1{((8*(eta2-1)*__nue - 4*eta2 + 7)*g - 2*eta2)/(8*c)};
    const _T c2{-((8*(eta2-1)*__nue - 4*eta2 + 1)*g + 2*eta2)/(16*c)};
    const _T c3{((24*(eta2-1)*__nue - 12*eta2 - 3)*g + 16*(1-eta2)*__nue + 10*eta2)/(16*c)};
    const _T c4{-((12*eta2+3)*g-10*eta2)/(16*c)};
    const _T c5{-((12*(eta2 - 1)*__nue + 15)*g+8*(1-eta2)*__nue - 2*eta2 - 8)/(8*c)};
    const _T c6{((60*eta2+45)*g - 54*eta2 - 16)/(16*c)};

    return tmech::eval(c1*II1 + c2*II2 + c3*II3 + c4*II4 + c5*II5 + c6*II6);
}

#endif // ESHELBY_TENSOR_SOLID_CYLINDRICAL_ELLIPSOID_MEAT_H
