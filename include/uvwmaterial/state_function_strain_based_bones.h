/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef STATE_FUNCTION_STRAIN_BASED_BONES_H
#define STATE_FUNCTION_STRAIN_BASED_BONES_H


template<typename T, std::size_t Dim, typename Container>
class von_mises_strain_state_function : public state_function_base<T, Dim, Container>
{
public:
    using value_type = T;

    von_mises_strain_state_function();

    von_mises_strain_state_function(solid_material_base<T, Dim, Container> & material);

    virtual ~von_mises_strain_state_function(){}

    inline void update() override;

    inline value_type value() const override;

    inline tmech::tensor<value_type, Dim, 2> derivative() const override;

};

template<typename T, std::size_t Dim, typename Container>
class vector_strain_state_function : public state_function_base<T, Dim, Container>
{
public:
    using value_type = T;

    vector_strain_state_function();

    template<typename Tensor>
    vector_strain_state_function(Tensor const& __m, Tensor const& __n);

    template<typename Tensor>
    vector_strain_state_function(solid_material_base<T, Dim, Container> & material, Tensor const& __m, Tensor const& __n);

    inline void update() override;

    inline value_type value() const override;

    inline tmech::tensor<value_type, Dim, 2> derivative() const override;

    template<typename Tensor>
    constexpr inline auto set_direction_vector(Tensor const& __m);

    template<typename Tensor>
    constexpr inline auto set_tangential_vector(Tensor const& __n);

private:
    tmech::tensor<value_type, Dim, 1> m;
    tmech::tensor<value_type, Dim, 1> n;
};


template<typename T, std::size_t Dim, typename Container>
class strain_energy_strain_state_function : public state_function_base<T, Dim, Container>
{
public:
    using value_type = T;

    strain_energy_strain_state_function();

    strain_energy_strain_state_function(value_type E);

    strain_energy_strain_state_function(solid_material_base<value_type, Dim, Container>  & material);

    strain_energy_strain_state_function(solid_material_base<value_type, Dim, Container>  & material, value_type E);

    constexpr inline void update() override;

    constexpr inline value_type value() const override;

    constexpr inline tmech::tensor<value_type, Dim, 2> derivative() const override;

    constexpr inline auto set_youngs_modulus(value_type const __E);

private:
    value_type E;
};


#endif // STATE_FUNCTION_STRAIN_BASED_BONES_H
