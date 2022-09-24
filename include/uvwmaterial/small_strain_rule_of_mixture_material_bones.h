/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef MATERIAL_SMALL_STRAIN_RULE_OF_MIXTURE_BONES_H
#define MATERIAL_SMALL_STRAIN_RULE_OF_MIXTURE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class small_strain_rule_of_mixture :
        public small_strain_material_base<_T, _Dim, _Container>,
        public rule_of_mixture_composite_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr small_strain_rule_of_mixture(){}

    virtual ~small_strain_rule_of_mixture(){}

    inline virtual void init() override {
        if(!small_strain_material_base<_T, _Dim, _Container>::_is_init || !rule_of_mixture_composite_base<_T, _Dim, _Container>::_is_init){
            rule_of_mixture_composite_base<_T, _Dim, _Container>::init();

            this->_C.fill(0);
            for(size_type i{0}; i<this->_materials.size(); ++i){
                this->_C += this->_volume_fractions[i]*tmech::dcontract(this->_materials[i]->tangent_tensor(), this->_strain_concentration_tensors[i]);
            }

            small_strain_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    inline virtual void reinit() override {
        rule_of_mixture_composite_base<_T, _Dim, _Container>::reinit();
        this->_C.fill(0);
        for(size_type i{0}; i<this->_materials.size(); ++i){
            this->_C += this->_volume_fractions[i]*tmech::dcontract(this->_materials[i]->tangent_tensor(), this->_strain_concentration_tensors[i]);
        }
    }

    inline virtual void update() override{
        const auto& strain{this->_strain};
        this->_C.fill(0);
        this->_stress.fill(0);
        for(size_type i{0}; i<this->_materials.size(); ++i){
            this->_materials[i]->strain_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], strain);
            this->_materials[i]->update();
            this->_C += this->_volume_fractions[i]*tmech::dcontract(this->_materials[i]->tangent_tensor(),this->_strain_concentration_tensors[i]);
            this->_stress += this->_volume_fractions[i]*this->_materials[i]->stress_tensor();
        }
    }

    inline virtual void update_stress() override{
        const auto& strain{this->_strain};
        this->_stress.fill(0);
        for(size_type i{0}; i<this->_materials.size(); ++i){
            this->_materials[i]->strain_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], strain);
            this->_materials[i]->update_stress();
            this->_stress += this->_volume_fractions[i]*this->_materials[i]->stress_tensor();
        }
    }

    inline virtual void update_tangent() override{
        const auto& strain{this->_strain};
        this->_C.fill(0);
        for(size_type i{0}; i<this->_materials.size(); ++i){
            this->_materials[i]->strain_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], strain);
            this->_materials[i]->update_tangent();
            this->_C += this->_volume_fractions[i]*tmech::dcontract(this->_materials[i]->tangent_tensor(), this->_strain_concentration_tensors[i]);
        }
    }

    inline virtual void update_strain()override{
        const auto& strain{this->_strain};
        for(size_type i{0}; i<this->_materials.size(); ++i){
            this->_materials[i]->strain_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], strain);
        }
    }

    inline virtual solid_material_base<_T, _Dim, _Container> * base_material() override {
        return this;
    }
};

#endif // MATERIAL_SMALL_STRAIN_RULE_OF_MIXTURE_BONES_H
