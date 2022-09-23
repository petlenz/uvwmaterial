#ifndef MEAN_FIELD_COMPOSITE_CONDUCTIVITY_BASE_BONES_H
#define MEAN_FIELD_COMPOSITE_CONDUCTIVITY_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_conductivity_base :
        public mean_field_composite_base<_T, _Dim, _Container>,
        public composite_material_conductivity_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    mean_field_composite_conductivity_base(){}

    virtual ~mean_field_composite_conductivity_base(){}

    virtual inline _T& volume_fraction(size_type __idx)override{
        return mean_field_composite_base<_T, _Dim, _Container>::volume_fraction(__idx);
    }

    virtual inline const _T& volume_fraction(size_type __idx) const override{
        return mean_field_composite_base<_T, _Dim, _Container>::volume_fraction(__idx);
    }

    virtual inline size_type number_of_materials() const override{
        return mean_field_composite_base<_T, _Dim, _Container>::number_of_materials();
    }

    virtual inline material_base<_T, _Dim, _Container> const* material(size_type __idx) const override{
        return mean_field_composite_base<_T, _Dim, _Container>::material(__idx);
    }

    virtual inline material_base<_T, _Dim, _Container>* material(size_type __idx)override{
        return mean_field_composite_base<_T, _Dim, _Container>::material(__idx);
    }

    virtual inline void init() override {
        if(!this->_is_init){
            this->check_data();
            composite_material_base<_T, _Dim, _Container>::init_materials();
            this->_gradient_concentration_tensors.resize(this->_inclusions.size()+1);
            dynamic_cast<mean_field_composite_conductivity_kernal_base<_T, _Dim, _Container>*>(this->_kernal)->determine_gradient_concentration_tensors();
            this->_is_init = true;
        }
    }

    virtual inline void reinit() override {
        dynamic_cast<mean_field_composite_conductivity_kernal_base<_T, _Dim, _Container>*>(this->_kernal)->determine_gradient_concentration_tensors();
    }
};

#endif // MEAN_FIELD_COMPOSITE_CONDUCTIVITY_BASE_BONES_H
