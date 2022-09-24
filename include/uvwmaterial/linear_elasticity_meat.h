/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef LINEAR_ELASTICITY_MEAT_H
#define LINEAR_ELASTICITY_MEAT_H

template <typename T, std::size_t Dim, typename Container>
linear_elasticity<T, Dim, Container>::linear_elasticity()
{}

template <typename T, std::size_t Dim, typename Container>
linear_elasticity<T, Dim, Container>::linear_elasticity(value_type const __E, value_type const __nu):
    small_strain_material_base<T ,Dim, Container>(__E, __nu)
{}

template <typename T, std::size_t Dim, typename Container>
linear_elasticity<T, Dim, Container>::linear_elasticity(std::vector<value_type> const& __parameter):
    small_strain_material_base<T ,Dim, Container>(__parameter)
{}

template <typename T, std::size_t Dim, typename Container>
linear_elasticity<T, Dim, Container>::linear_elasticity(std::initializer_list<value_type> const& __parameter):
    small_strain_material_base<T ,Dim, Container>(__parameter)
{}

template <typename T, std::size_t Dim, typename Container>
template<typename ...Parameter>
linear_elasticity<T, Dim, Container>::linear_elasticity(Parameter ... __parameter):
    small_strain_material_base<T ,Dim, Container>(__parameter...)
{}

template <typename T, std::size_t Dim, typename Container>
inline void linear_elasticity<T, Dim, Container>::init(){
    if(!this->_is_init){
        if(this->_parameter.empty()){
            throw std::runtime_error("linear_elasticity::init() no matching number of properties");
        }

        constexpr value_type fac{static_cast<value_type>(Dim)};
        const auto I{tmech::eye<value_type, Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
        const auto IIvol{tmech::otimes(I, I)/fac};
        const auto IIdev{IIsym - IIvol};

        const auto [_E, _nu]{get_parameter_base()};
        const value_type _K{K_E_nu(_E,_nu)};
        const value_type _mu{mu_E_nu(_E,_nu)};
        this->_C = 3*_K*IIvol + 2*_mu*IIdev;
        //(fac/fac) --> IIvol --> II

        this->_is_init = true;
    }
}

template <typename T, std::size_t Dim, typename Container>
inline void linear_elasticity<T, Dim, Container>::update(){
    const auto I{tmech::eye<value_type, Dim, 2>()};
    const auto [_E, _nu]{get_parameter_base()};
    const value_type _lambda{uvwmat::lambda_E_nu(_E,_nu)};
    const value_type _mu{uvwmat::mu_E_nu(_E,_nu)};
    this->_stress = _lambda*I*tmech::trace(this->strain_tensor()) + 2.0*_mu*this->strain_tensor();
    //this->_stress = tmech::dcontract(this->_C, this->strain_tensor());
}

template <typename T, std::size_t Dim, typename Container>
inline void linear_elasticity<T, Dim, Container>::update_stress(){
    this->_stress = tmech::dcontract(this->_C, this->strain_tensor());
}

template <typename T, std::size_t Dim, typename Container>
inline void linear_elasticity<T, Dim, Container>::update_tangent(){
    //nothing to do
}

#endif // LINEAR_ELASTICITY_MEAT_H
