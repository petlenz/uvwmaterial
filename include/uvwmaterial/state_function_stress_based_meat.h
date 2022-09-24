/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef STATE_FUNCTION_STRESS_BASED_MEAT_H
#define STATE_FUNCTION_STRESS_BASED_MEAT_H


template<typename T, std::size_t Dim, typename Container>
von_mises_state_function<T, Dim, Container>::von_mises_state_function():
    state_function_base<T, Dim, Container>()
{}

template<typename T, std::size_t Dim, typename Container>
von_mises_state_function<T, Dim, Container>::von_mises_state_function(solid_material_base<T, Dim, Container> & material):
    state_function_base<T, Dim, Container>(material)
{}

template<typename T, std::size_t Dim, typename Container>
inline void von_mises_state_function<T, Dim, Container>::update(){
    assert(this->_material != nullptr);
    const auto& material{*this->_material};
    auto sig_dev{tmech::dev(material.stress_tensor())};
    this->_state_value = std::sqrt(1.5*tmech::dcontract(sig_dev,sig_dev));
}

template<typename T, std::size_t Dim, typename Container>
inline typename von_mises_state_function<T, Dim, Container>::value_type von_mises_state_function<T, Dim, Container>::value() const {
    return this->_state_value;
}

template<typename T, std::size_t Dim, typename Container>
inline tmech::tensor<typename von_mises_state_function<T, Dim, Container>::value_type, Dim, 2> von_mises_state_function<T, Dim, Container>::derivative() const {
    assert(this->_material != nullptr);
    const auto& material{*this->_material};
    return static_cast<tmech::tensor<T, Dim, 2>>(1.5*dev(material.stress_tensor())/(this->_state_value));
}


#endif // STATE_FUNCTION_STRESS_BASED_MEAT_H
