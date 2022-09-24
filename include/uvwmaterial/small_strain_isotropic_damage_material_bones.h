/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SMALL_STRAIN_ISOTROPIC_DAMAGE_MATERIAL_BONES_H
#define SMALL_STRAIN_ISOTROPIC_DAMAGE_MATERIAL_BONES_H

template <typename T, std::size_t Dim, typename Container>
class small_strain_isotropic_damage :
        public history_material_base<T>,
        public small_strain_material_base<T, Dim, Container>
{
public:
    using value_type = T;
    using size_type = std::size_t;
    using tensor4 = tmech::tensor<value_type, Dim, 4>;

    small_strain_isotropic_damage();

    small_strain_isotropic_damage(small_strain_material_base<T, Dim, Container> * __material,
                                  propagation_law_base<value_type> * __damage_func,
                                  yield_function_isotropic_damage_strain_based<T, Dim, Container> * __yield_function);

    virtual ~small_strain_isotropic_damage(){}

    inline virtual void init()override;

    inline virtual void reinit()override{
        _material->reinit();
        this->_C = _material->tangent_tensor();
    }

    inline virtual void update()override;

    inline virtual void update_stress()override;

    inline virtual void update_tangent()override;

    constexpr inline auto base_material(small_strain_material_base<value_type, Dim, Container> * __material);

    constexpr inline auto propagation_law(propagation_law_base<value_type> * __propagation_law);

    constexpr inline auto yield_function(yield_function_isotropic_damage_strain_based<value_type, Dim, Container> * __yield_function);

    virtual inline solid_material_base<T, Dim, Container>* base_material();

private:
    constexpr inline auto check_data(){
        if(_material == nullptr){
            throw std::runtime_error("small_strain_isotropic_damage::init(): material not set");
        }
        if(_yield_function == nullptr){
            throw std::runtime_error("small_strain_isotropic_damage::init(): yield function not set");
        }
        if(_damage_func == nullptr){
            throw std::runtime_error("small_strain_isotropic_damage::init(): propagation law not set");
        }
    }

    small_strain_material_base<value_type, Dim, Container> * _material;
    propagation_law_base<value_type> * _damage_func;
    yield_function_isotropic_damage_strain_based<value_type, Dim, Container> * _yield_function;
};

#endif // SMALL_STRAIN_ISOTROPIC_DAMAGE_MATERIAL_BONES_H
