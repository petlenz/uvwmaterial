/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SMALL_STRAIN_ISOTROPIC_DAMAGE_MATERIAL_MEAT_H
#define SMALL_STRAIN_ISOTROPIC_DAMAGE_MATERIAL_MEAT_H

template <typename T, std::size_t Dim, typename Container>
small_strain_isotropic_damage<T, Dim, Container>::small_strain_isotropic_damage():
    history_material_base<T>(2),
    small_strain_material_base<T, Dim, Container>(),
    _material(nullptr),
    _yield_function(nullptr)
{}

template <typename T, std::size_t Dim, typename Container>
small_strain_isotropic_damage<T, Dim, Container>::small_strain_isotropic_damage(small_strain_material_base<T, Dim, Container> * __material,
                                                                                propagation_law_base<value_type> * __damage_func,
                                                                                yield_function_isotropic_damage_strain_based<T, Dim, Container> * __yield_function):
    history_material_base<T>(2),
    small_strain_material_base<T, Dim, Container>(),
    _material(__material),
    _damage_func(__damage_func),
    _yield_function(__yield_function)
{}

template <typename T, std::size_t Dim, typename Container>
inline void small_strain_isotropic_damage<T, Dim, Container>::init(){
    if(!this->_is_init){
        check_data();
        _material->init();
        this->_C = _material->tangent_tensor();
        this->_is_init = true;
    }
}

template <typename T, std::size_t Dim, typename Container>
inline void small_strain_isotropic_damage<T, Dim, Container>::update(){
    //easy handling
    auto& sig{this->_stress};
    auto& C{this->_C};
    auto& history{this->_history};

    //update base material
    //compute trail stress
    _material->strain_tensor() = this->_strain;

    //set history
    auto* material_his{dynamic_cast<history_material_base<T>*>(_material)};
    if(material_his){
        material_his->set_history(&history[2], &history[history.size()]);
    }

    //set temperature

    //set eigenstrain....

    //update stress and strains of basis material
    _material->update();

    //get history
    if(material_his){
        for(size_type i{0}; i<history.size()-2; ++i){
            history[i+2] = material_his->history()[i];
        }
    }

    //update equivalent scalar
    _yield_function->update_equivalent_scalar();

    //solve for new damage multiplier
    const auto his_new{_yield_function->solve(history[1])};

    //check if yielding
    if(_yield_function->yielding(history[1])){
        //Damage
        //std::cout<<" "<<_yield_function->state()<<" "<<_damage_func->value(his_new)<<" "<<_damage_func->derivative(his_new)<<std::endl;
        //std::cout<<his_new<<std::endl;
        _yield_function->update_strain_dependent_parts();
        const auto D{_damage_func->value(his_new)};
        const auto dD{_damage_func->derivative(his_new)};
        const auto& dEps_eq{_yield_function->derivative()};
        C = (1-D)*_material->tangent_tensor() - dD*otimes(_material->stress_tensor(), dEps_eq);
        sig = (1-D)*_material->stress_tensor();
        history[0] = D;
    }else{
        //No damage
        const auto D{history[0]};
        C = (1-D)*_material->tangent_tensor();
        sig = (1-D)*_material->stress_tensor();
    }
    history[1] = his_new;
}

template <typename T, std::size_t Dim, typename Container>
inline void small_strain_isotropic_damage<T, Dim, Container>::update_stress(){
    throw std::runtime_error("not implemented yet");
}

template <typename T, std::size_t Dim, typename Container>
inline void small_strain_isotropic_damage<T, Dim, Container>::update_tangent(){
    throw std::runtime_error("not implemented yet");
}

template <typename T, std::size_t Dim, typename Container>
constexpr inline auto small_strain_isotropic_damage<T, Dim, Container>::base_material(small_strain_material_base<value_type, Dim, Container> * __material){
    _material = __material;
}

template <typename T, std::size_t Dim, typename Container>
constexpr inline auto small_strain_isotropic_damage<T, Dim, Container>::propagation_law(propagation_law_base<value_type> * __propagation_law){
    _damage_func = __propagation_law;
}

template <typename T, std::size_t Dim, typename Container>
constexpr inline auto small_strain_isotropic_damage<T, Dim, Container>::yield_function(yield_function_isotropic_damage_strain_based<value_type, Dim, Container> * __yield_function){
    _yield_function = __yield_function;
}

template <typename T, std::size_t Dim, typename Container>
inline solid_material_base<T, Dim, Container>* small_strain_isotropic_damage<T, Dim, Container>::base_material(){
    return _material->base_material();
}

#endif // SMALL_STRAIN_ISOTROPIC_DAMAGE_MATERIAL_MEAT_H
