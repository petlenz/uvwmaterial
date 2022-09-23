#ifndef MEAN_FIELD_COMPOSITE_CONDUCTIVITY_MATERIAL_BONES_H
#define MEAN_FIELD_COMPOSITE_CONDUCTIVITY_MATERIAL_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_conductivity_material :
        public mean_field_composite_conductivity_base<_T, _Dim, _Container>,
        public conductivity_material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    mean_field_composite_conductivity_material() {}

    virtual ~mean_field_composite_conductivity_material() {}

    virtual inline void init() override {
        if(!mean_field_composite_conductivity_base<_T, _Dim, _Container>::_is_init || !conductivity_material_base<_T, _Dim, _Container>::_is_init){
            this->_matrix_material->init();
            for(auto inclusion : this->_inclusions){
                inclusion->init();
            }

            //new
            _materials.reserve(this->_inclusions.size()+1);
            for(size_type i{0}; i<this->_inclusions.size(); ++i){
                _materials.push_back(make_conductivity_material(this->_inclusions[i]->material()));
            }
            _materials.push_back(make_conductivity_material(this->_matrix_material));

            mean_field_composite_conductivity_base<_T, _Dim, _Container>::init();

            for(size_type i{0}; i<this->_materials.size(); ++i){
                const auto ci{this->volume_fraction(i)};
                this->_conductivity += ci*_materials[i]->conductivity_tensor()*this->_gradient_concentration_tensors[i];
            }

            conductivity_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    virtual inline void reinit() override {
        this->_matrix_material->reinit();
        for(auto inclusion : this->_inclusions){
            inclusion->reinit();
        }
        mean_field_composite_conductivity_base<_T, _Dim, _Container>::reinit();
        for(size_type i{0}; i<this->_materials.size(); ++i){
            const auto ci{this->volume_fraction(i)};
            this->_conductivity += ci*_materials[i]->conductivity_tensor()*this->_gradient_concentration_tensors[i];
        }
    }

    virtual inline void update_gradient()override{
        const auto& gradient_macro{this->gradient_tensor()};
        for(size_type i{0}; i<this->_materials.size(); ++i){
            _materials[i]->gradient_tensor() = this->_gradient_concentration_tensors[i]*gradient_macro;
        }
    }

    virtual inline void update() override {
        update_gradient();
        this->_flux.fill(0);
        this->_conductivity.fill(0);
        for(size_type i{0}; i<this->_materials.size(); ++i){
            _materials[i]->update();
            const auto ci{this->volume_fraction(i)};
            this->_flux += ci*_materials[i]->flux_tensor();
            this->_conductivity += ci*_materials[i]->conductivity_tensor()*this->_gradient_concentration_tensors[i];
        }
    }

    virtual inline void update_flux() override {
        update_gradient();
        this->_flux.fill(0);
        for(size_type i{0}; i<this->_materials.size(); ++i){
            _materials[i]->update_flux();
            const auto ci{this->volume_fraction(i)};
            this->_flux += ci*_materials[i]->flux_tensor();
        }
    }

    virtual inline void update_conductivity() override {
        update_gradient();
        this->_conductivity.fill(0);
        for(size_type i{0}; i<this->_materials.size(); ++i){
            _materials[i]->update_conductivity();
            const auto ci{this->volume_fraction(i)};
            this->_conductivity += ci*_materials[i]->conductivity_tensor()*this->_gradient_concentration_tensors[i];
        }
    }

    virtual inline conductivity_material_base<_T, _Dim, _Container>* base_material() override {
        return this;
    }

protected:
    std::vector<conductivity_material_base<_T, _Dim, _Container>*> _materials;
};


#endif // MEAN_FIELD_COMPOSITE_CONDUCTIVITY_MATERIAL_BONES_H
