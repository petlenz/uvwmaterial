#ifndef SMALL_STRAIN_MEAN_FIELD_COMPOSITE_BONES_H
#define SMALL_STRAIN_MEAN_FIELD_COMPOSITE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class small_strain_mean_field_composite :
        public mean_field_composite_solid_base<_T, _Dim, _Container>,
        public small_strain_material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    small_strain_mean_field_composite() {}

    virtual ~small_strain_mean_field_composite() {}

    inline virtual void init() override {
        if(!small_strain_material_base<_T, _Dim, _Container>::_is_init || !mean_field_composite_solid_base<_T, _Dim, _Container>::_is_init){
            this->_matrix_material->init();
            for(auto inclusion : this->_inclusions){
                inclusion->init();
            }

            mean_field_composite_solid_base<_T, _Dim, _Container>::init();
            auto matrix_material = make_solid_material(this->_matrix_material);
            this->_C = this->_cm*tmech::dcontract(matrix_material->tangent_tensor(),this->_strain_concentration_tensors.back());

            for(size_type i{0}; i<this->_inclusions.size(); ++i){
                const auto ci{this->_inclusions[i]->volume_fraction()};
                this->_C += ci*tmech::dcontract(make_solid_material(this->_inclusions[i]->material())->tangent_tensor(), this->_strain_concentration_tensors[i]);
            }

            //new
            _matrix_material_small_strain = dynamic_cast<small_strain_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(size_type i{0}; i<this->_inclusions.size(); ++i){
                _inclusions_material.push_back(dynamic_cast<small_strain_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material()));
            }


            small_strain_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    inline virtual void reinit() override {
        this->_matrix_material->reinit();
        for(auto inclusion : this->_inclusions){
            inclusion->reinit();
        }
        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();
        this->_C = this->_cm*tmech::dcontract(_matrix_material_small_strain->tangent_tensor(), this->_strain_concentration_tensors.back());
        for(size_type i{0}; i<this->_inclusions.size(); ++i){
            this->_C += this->_inclusions[i]->volume_fraction()*tmech::dcontract(_inclusions_material[i]->tangent_tensor(), this->_strain_concentration_tensors[i]);
        }
    }

    inline virtual void update_strain()override{
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_strain(): class is not initialized");
        }

        const auto& eps_macro{this->strain_tensor()};
        _matrix_material_small_strain->strain_tensor() = tmech::dcontract(this->_strain_concentration_tensors.back(), eps_macro);
        for(size_type i{0}; i<this->_inclusions.size(); ++i){
            _inclusions_material[i]->strain_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], eps_macro);
        }
    }

    inline virtual void update() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update(): is not initialized");
        }

        update_strain();

        //auto matrix_material{dynamic_cast<small_strain_material_base<_T, _Dim, _Container>*>(this->_matrix_material)};
        _matrix_material_small_strain->update();

        this->_stress = this->_cm*_matrix_material_small_strain->stress_tensor();
        this->_C = this->_cm*tmech::dcontract(_matrix_material_small_strain->tangent_tensor(), this->_strain_concentration_tensors.back());
        for(size_type i{0}; i<this->_inclusions.size(); ++i){
            auto small_strain_material{_inclusions_material[i]};
            small_strain_material->update();
            this->_stress += this->_inclusions[i]->volume_fraction()*small_strain_material->stress_tensor();
            this->_C += this->_inclusions[i]->volume_fraction()*tmech::dcontract(small_strain_material->tangent_tensor(), this->_strain_concentration_tensors[i]);
        }
    }

    inline virtual void update_stress() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_stress(): is not initialized");
        }

        update_strain();
        _matrix_material_small_strain->update_stress();
        this->_stress = this->_cm*_matrix_material_small_strain->stress_tensor();

        for(size_type i{0}; i<this->_inclusions.size(); ++i){
            auto small_strain_material{_inclusions_material[i]};
            small_strain_material->update_stress();
            this->_stress += this->_inclusions[i]->volume_fraction()*small_strain_material->stress_tensor();
        }
    }

    inline virtual void update_tangent() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_tangent(): is not initialized");
        }

        update_strain();
        _matrix_material_small_strain->update_tangent();
        this->_C = this->_cm*dcontract(_matrix_material_small_strain->tangent_tensor(), this->_strain_concentration_tensors.back());

        for(size_type i{0}; i<this->_inclusions.size(); ++i){
            auto small_strain_material{_inclusions_material[i]};
            small_strain_material->update_tangent();
            this->_C += this->_inclusions[i]->volume_fraction()*dcontract(small_strain_material->tangent_tensor(), this->_strain_concentration_tensors[i]);
        }
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override {
        return this;
    }


private:
    small_strain_material_base<_T, _Dim, _Container>* _matrix_material_small_strain;
    std::vector<small_strain_material_base<_T, _Dim, _Container>*> _inclusions_material;
};

#endif // SMALL_STRAIN_MEAN_FIELD_COMPOSITE_BONES_H
