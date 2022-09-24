/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef PLASTIC_MATERIAL_ISOTROPIC_YIELDING_BASE_BONES_H
#define PLASTIC_MATERIAL_ISOTROPIC_YIELDING_BASE_BONES_H

template <typename T, std::size_t Dim, typename Container>
class isotropic_yield_function_base
{
public:
    using value_type = T;

    isotropic_yield_function_base():
        history_n(0),
        scalar_eq(0),
        critical_val(0),
        parameter(),
        base_material(nullptr)
    {}

    isotropic_yield_function_base(material_base<T, Dim, Container> * __base_material):
        history_n(0),
        scalar_eq(0),
        critical_val(0),
        parameter(),
        base_material(__base_material)
    {}

    template<typename ...Parameter>
    isotropic_yield_function_base(material_base<T, Dim, Container> * __base_material, value_type __critical_val, Parameter && ... __parameter):
        history_n(0),
        scalar_eq(0),
        critical_val(__critical_val),
        parameter(),
        base_material(__base_material)
    {
        set_parameter(__parameter...);
    }

    virtual ~isotropic_yield_function_base(){}

    constexpr inline auto set_base_material(material_base<T, Dim, Container> * __material_base){
        base_material = __material_base;
    }

    constexpr inline auto set_state_function(state_function_base<T, Dim, Container> * __state_function){
        equivalent_func = __state_function;
    }

    constexpr inline auto set_critical_value(value_type const __critical_val){
        critical_val = __critical_val;
    }

    template<typename ...Parameter>
    constexpr inline auto set_parameter(Parameter && ... __parameter){
        set_parameter_in_container<std::vector<value_type>>::set_parameter(parameter, __parameter...);
    }

    virtual constexpr inline bool yielding(value_type const& __history) const = 0;

    virtual constexpr inline value_type solve(value_type const& __history) = 0;

    constexpr inline auto get_state()const{
        return scalar_eq;
    }

    constexpr inline auto get_critical()const{
        return critical_val;
    }

    constexpr inline void update_equivalent_scalar(){
        assert(equivalent_func != nullptr);
        equivalent_func->update();
        scalar_eq = equivalent_func->value();
    }

protected:
    value_type history_n;
    value_type scalar_eq;
    value_type critical_val;
    std::vector<value_type> parameter;
    state_function_base<T, Dim, Container> * equivalent_func;
    material_base<T, Dim, Container> * base_material;
};



#endif // PLASTIC_MATERIAL_ISOTROPIC_YIELDING_BASE_BONES_H
