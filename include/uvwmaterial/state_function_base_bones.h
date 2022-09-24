/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef STATE_FUNCTION_BASE_BONES_H
#define STATE_FUNCTION_BASE_BONES_H

template <typename T, std::size_t Dim, typename Container>
class state_function_base
{
public:
    using value_type = T;

    state_function_base();

    state_function_base(solid_material_base<value_type, Dim, Container> * _material);

    virtual ~state_function_base(){}

    virtual inline void update() = 0;

    virtual inline value_type value() const = 0;

    virtual inline tmech::tensor<value_type, Dim, 2> derivative() const = 0;

    constexpr inline auto set_base_material(solid_material_base<value_type, Dim, Container> * __material);

protected:
    value_type _state_value;
    solid_material_base<value_type, Dim, Container> * _material;
};


//template <std::size_t Dim, typename Container>
//class state_function_base<float, Dim, Container>
//{
//public:
//    using value_type = float;

//    state_function_base();

//    state_function_base(material_base<value_type, Dim, Container> * _material);

//    virtual inline void update() = 0;

//    virtual inline value_type value() const = 0;

//    virtual inline tmech::tensor<value_type, Dim, 2> derivative() const = 0;

//    constexpr inline auto set_base_material(material_base<value_type, Dim, Container> * __material);

//protected:
//    value_type state_value;
//    material_base<value_type, Dim, Container> * material;
//};

#endif // STATE_FUNCTION_BASE_BONES_H
