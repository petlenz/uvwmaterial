/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef STATE_FUNCTION_STRAIN_BASED_MEAT_H
#define STATE_FUNCTION_STRAIN_BASED_MEAT_H


template<typename T, std::size_t Dim, typename Container>
von_mises_strain_state_function<T, Dim, Container>::von_mises_strain_state_function():
    state_function_base<T, Dim, Container>()
{}

template<typename T, std::size_t Dim, typename Container>
von_mises_strain_state_function<T, Dim, Container>::von_mises_strain_state_function(solid_material_base<T, Dim, Container> & material):
    state_function_base<T, Dim, Container>(material)
{}

template<typename T, std::size_t Dim, typename Container>
inline void von_mises_strain_state_function<T, Dim, Container>::update(){
    assert(this->_material != nullptr);
    const auto& material{*make_small_strain_material(this->_material)};
    this->_state_value = std::sqrt(1.5*tmech::dcontract(tmech::dev(material.strain_tensor()),tmech::dev(material.strain_tensor())));
}

template<typename T, std::size_t Dim, typename Container>
inline typename von_mises_strain_state_function<T, Dim, Container>::value_type von_mises_strain_state_function<T, Dim, Container>::value() const {
    return this->_state_value;
}

template<typename T, std::size_t Dim, typename Container>
inline tmech::tensor<typename von_mises_strain_state_function<T, Dim, Container>::value_type, Dim, 2> von_mises_strain_state_function<T, Dim, Container>::derivative() const {
    assert(this->_material != nullptr);
    const auto& material{*make_small_strain_material(this->_material)};
    return static_cast<tmech::tensor<T, Dim, 2>>(1.5*tmech::dev(material.strain_tensor())/(this->_state_value));
}

template<typename T, std::size_t Dim, typename Container>
vector_strain_state_function<T, Dim, Container>::vector_strain_state_function():
    state_function_base<T, Dim, Container>(),
    m(),
    n()
{}

template<typename T, std::size_t Dim, typename Container>
template<typename Tensor>
vector_strain_state_function<T, Dim, Container>::vector_strain_state_function(Tensor const& __m, Tensor const& __n):
    state_function_base<T, Dim, Container>(),
    m(__m),
    n(__n)
{}

template<typename T, std::size_t Dim, typename Container>
template<typename Tensor>
vector_strain_state_function<T, Dim, Container>::vector_strain_state_function(solid_material_base<T, Dim, Container> & material, Tensor const& __m, Tensor const& __n):
    state_function_base<T, Dim, Container>(material),
    m(__m),
    n(__n)
{}

template<typename T, std::size_t Dim, typename Container>
inline void vector_strain_state_function<T, Dim, Container>::update() {
    assert(this->_material != nullptr);
    const auto& material{*make_small_strain_material(this->_material)};
    this->_state_value = tmech::dot(m, material.strain_tensor()*n);
    //this->_state_value = tmech::dot(m, material.stress_tensor()*n);
}

template<typename T, std::size_t Dim, typename Container>
inline typename vector_strain_state_function<T, Dim, Container>::value_type vector_strain_state_function<T, Dim, Container>::value() const {
    return std::abs(this->_state_value);
}

template<typename T, std::size_t Dim, typename Container>
inline tmech::tensor<typename vector_strain_state_function<T, Dim, Container>::value_type, Dim, 2> vector_strain_state_function<T, Dim, Container>::derivative() const {
    assert(this->_material != nullptr);
    //const auto& material{*make_small_strain_material(this->_material)};
    return static_cast<tmech::tensor<T, Dim, 2>>(tmech::sym(otimes(m,n))*tmech::sign(this->_state_value));
    //return static_cast<tmech::tensor<T, Dim, 2>>(tmech::dcontract(tmech::otimes(m,n), material.tangent_tensor())*tmech::sign(this->_state_value));
    //return static_cast<tmech::tensor<T, Dim, 2>>(0.5*(otimes(m,n) + otimes(n,m))*this->_state_value/std::abs(this->_state_value));
}

template<typename T, std::size_t Dim, typename Container>
template<typename Tensor>
constexpr inline auto vector_strain_state_function<T, Dim, Container>::set_direction_vector(Tensor const& __m){
    m = __m;
}

template<typename T, std::size_t Dim, typename Container>
template<typename Tensor>
constexpr inline auto vector_strain_state_function<T, Dim, Container>::set_tangential_vector(Tensor const& __n){
    n = __n;
}






template<typename T, std::size_t Dim, typename Container>
strain_energy_strain_state_function<T, Dim, Container>::strain_energy_strain_state_function():
    state_function_base<T, Dim, Container>(),
    E(0)
{}

template<typename T, std::size_t Dim, typename Container>
strain_energy_strain_state_function<T, Dim, Container>::strain_energy_strain_state_function(solid_material_base<value_type, Dim, Container>  & material):
    state_function_base<T, Dim, Container>(material),
    E(0)
{}

template<typename T, std::size_t Dim, typename Container>
strain_energy_strain_state_function<T, Dim, Container>::strain_energy_strain_state_function(value_type E):
    state_function_base<T, Dim, Container>(),
    E(E)
{}

template<typename T, std::size_t Dim, typename Container>
strain_energy_strain_state_function<T, Dim, Container>::strain_energy_strain_state_function(solid_material_base<value_type, Dim, Container> & material, value_type E):
    state_function_base<T, Dim, Container>(material),
    E(E)
{}

template<typename T, std::size_t Dim, typename Container>
constexpr inline void strain_energy_strain_state_function<T, Dim, Container>::update(){
    assert(this->_material != nullptr);
    const auto& material{*make_small_strain_material(this->_material)};
    this->_state_value = std::sqrt(tmech::dcontract(material.strain_tensor(), tmech::dcontract(material.tangent_tensor(),material.strain_tensor()))/E);
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline typename strain_energy_strain_state_function<T, Dim, Container>::value_type strain_energy_strain_state_function<T, Dim, Container>::value() const {
    //std::cout<<"von Mises state "<<this->_state_value<<std::endl;
    return this->_state_value;
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline tmech::tensor<typename strain_energy_strain_state_function<T, Dim, Container>::value_type, Dim, 2> strain_energy_strain_state_function<T, Dim, Container>::derivative() const {
    assert(this->_material != nullptr);
    const auto& material{*make_small_strain_material(this->_material)};
    return static_cast<tmech::tensor<T, Dim, 2>>(tmech::dcontract(material.tangent_tensor(),material.strain_tensor())/(E*this->_state_value));
}

template<typename T, std::size_t Dim, typename Container>
constexpr inline auto strain_energy_strain_state_function<T, Dim, Container>::set_youngs_modulus(value_type const __E){
    E = __E;
}


#endif // STATE_FUNCTION_STRAIN_BASED_MEAT_H
