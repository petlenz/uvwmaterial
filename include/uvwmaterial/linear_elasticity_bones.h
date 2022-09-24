/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef LINEAR_ELASTICITY_BONES_H
#define LINEAR_ELASTICITY_BONES_H

template <typename T, std::size_t Dim, typename Container>
class linear_elasticity : public small_strain_material_base<T, Dim, Container>
{
public:
    using value_type = T;
    using size_type = std::size_t;
    using tensor4 = tmech::tensor<value_type, Dim, 4>;

    linear_elasticity();

    linear_elasticity(value_type const __E, value_type const __nu);

    linear_elasticity(std::vector<value_type> const& __parameter);

    linear_elasticity(std::initializer_list<value_type> const& __parameter);

    template<typename ...Parameter>
    linear_elasticity(Parameter ... __parameter);

    virtual ~linear_elasticity(){}

    inline virtual void init() override;

    inline virtual void reinit() override{
        this->_is_init = false;
        this->init();
    }

    inline virtual void update() override;

    inline virtual void update_stress() override;

    inline virtual void update_tangent() override;

    inline virtual solid_material_base<T, Dim, Container>* base_material() override{
        return this;
    }

private:
    constexpr inline auto get_parameter_base()const{
        constexpr auto is_container_map{detail::is_map_container<Container>::value};
        constexpr auto is_container_vector_any{detail::is_vector_any_container<Container>::value};
        constexpr auto is_container_vector_value_type{detail::is_vector_value_type_container<Container>::value};
        static_assert (is_container_map || is_container_vector_any || is_container_vector_value_type, "NO MATCHING CONTAINER");

        if constexpr (is_container_map){
            return std::make_pair(this->_parameter["E"], this->_parameter["nu"]);
        }else if constexpr (is_container_vector_any){
            return std::make_pair(std::any_cast<value_type>(this->_parameter[0]), std::any_cast<value_type>(this->parameter[1]));
        }else if constexpr (is_container_vector_value_type) {
            return std::make_pair(this->_parameter[0], this->_parameter[1]);
        }
    }
};

#endif // LINEAR_ELASTICITY_BONES_H
