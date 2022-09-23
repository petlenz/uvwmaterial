#ifndef INCREMENTAL_incremental_linear_elasticity_MEAT_H
#define INCREMENTAL_incremental_linear_elasticity_MEAT_H

template <typename T, std::size_t Dim, typename Container>
incremental_linear_elasticity<T, Dim, Container>::incremental_linear_elasticity()
{}

template <typename T, std::size_t Dim, typename Container>
incremental_linear_elasticity<T, Dim, Container>::incremental_linear_elasticity(value_type const __E, value_type const __nu):
    small_strain_material_base<T ,Dim, Container>(__E, __nu)
{}

template <typename T, std::size_t Dim, typename Container>
incremental_linear_elasticity<T, Dim, Container>::incremental_linear_elasticity(std::vector<value_type> const& __parameter):
    small_strain_material_base<T ,Dim, Container>(__parameter)
{}

template <typename T, std::size_t Dim, typename Container>
incremental_linear_elasticity<T, Dim, Container>::incremental_linear_elasticity(std::initializer_list<value_type> const& __parameter):
    small_strain_material_base<T ,Dim, Container>(__parameter)
{}

template <typename T, std::size_t Dim, typename Container>
template<typename ...Parameter>
incremental_linear_elasticity<T, Dim, Container>::incremental_linear_elasticity(Parameter ... __parameter):
    small_strain_material_base<T ,Dim, Container>(__parameter...)
{}

template <typename T, std::size_t Dim, typename Container>
inline void incremental_linear_elasticity<T, Dim, Container>::init(){
    if(!this->_is_init){
        if(this->_parameter.empty()){
            throw std::runtime_error("incremental_linear_elasticity::init() no matching number of properties");
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
inline void incremental_linear_elasticity<T, Dim, Container>::reinit(){
    this->_is_init = false;
    this->init();
}

template <typename T, std::size_t Dim, typename Container>
inline void incremental_linear_elasticity<T, Dim, Container>::update(){
    const auto I{tmech::eye<value_type, Dim, 2>()};
    const auto [_E, _nu]{get_parameter_base()};
    const value_type _lambda{uvwmat::lambda_E_nu(_E,_nu)};
    const value_type _mu{uvwmat::mu_E_nu(_E,_nu)};
    this->_stress += _lambda*I*tmech::trace(this->dstrain_tensor()) + 2.0*_mu*this->dstrain_tensor();
}

template <typename T, std::size_t Dim, typename Container>
inline void incremental_linear_elasticity<T, Dim, Container>::update_stress(){
    this->_stress += dcontract(this->_C, this->dstrain_tensor());
}

template <typename T, std::size_t Dim, typename Container>
inline void incremental_linear_elasticity<T, Dim, Container>::update_tangent(){
    //nothing to do
}

#endif // INCREMENTAL_LINEAR_ELASTICITY_MEAT_H
