#ifndef MEAN_FIELD_COMPOSITE_CONDUCTIVITY_KERNAL_BASE_BONES_H
#define MEAN_FIELD_COMPOSITE_CONDUCTIVITY_KERNAL_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_conductivity_base;

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_conductivity_kernal_base :
        public mean_field_composite_kernal_base<_T, _Dim, _Container>
{
public:
    mean_field_composite_conductivity_kernal_base(){}

    virtual ~mean_field_composite_conductivity_kernal_base(){}

    inline virtual void determine_gradient_concentration_tensors() = 0;
};


//dilute
template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_conductivity_dilute_kernal : public mean_field_composite_conductivity_kernal_base<_T, _Dim, _Container>
{
public:
    using size_type  = std::size_t;
    using value_type = _T;

    mean_field_composite_conductivity_dilute_kernal(){}

    virtual ~mean_field_composite_conductivity_dilute_kernal() {}


    inline virtual void determine_gradient_concentration_tensors() override {
        auto& A{dynamic_cast<composite_material_conductivity_base<_T, _Dim, _Container>*>(this->_material)->gradient_concentration_tensors()};
        const auto& inclusions{this->_material->inclusions()};
        const auto* matrix_material{make_conductivity_material(this->_material->matrix_material())};

        const tmech::eye<value_type, _Dim, 2> I;
        const auto& Cm{matrix_material->conductivity_tensor()};
        A.back() = I;

        const tmech::tensor<value_type, _Dim, 2> Cm_inv{tmech::inv(Cm)};
        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto& Ci{make_conductivity_material(inclusions[i]->material())->conductivity_tensor()};
            const auto& S{dynamic_cast<eshelby_tensor_conductivity_base<_T, _Dim>*>(inclusions[i]->eshelby_tensor())->tensor()};
            A[i] = tmech::inv(I + S*Cm_inv*(Ci - Cm));
        }
    }
};

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_conductivity_mori_tanaka_kernal : public mean_field_composite_conductivity_kernal_base<_T, _Dim, _Container>
{
public:
    using size_type  = std::size_t;
    using value_type = _T;

    mean_field_composite_conductivity_mori_tanaka_kernal() {}

    inline virtual void determine_gradient_concentration_tensors() override {
        auto& A{dynamic_cast<composite_material_conductivity_base<_T, _Dim, _Container>*>(this->_material)->gradient_concentration_tensors()};
        const auto& inclusions{this->_material->inclusions()};
        const auto* matrix_material{make_conductivity_material(this->_material->matrix_material())};
        std::vector<tmech::tensor<_T, _Dim, 2>> TT(inclusions.size());

        tmech::tensor<value_type, _Dim, 2> TT_sum;
        const tmech::eye<value_type, _Dim, 2> I;

        //A_i = A_dil_i * inv(sum_r c_r A_dil_r)
        //r in [0, n+1]; i in [0, n+1]
        //A_dil_n+1 = cm*I;

        TT_sum = I*this->_material->matrix_volume_fraction();

        const auto& Cm{matrix_material->conductivity_tensor()};
        const tmech::tensor<value_type, _Dim, 2> Cm_inv{tmech::inv(Cm)};
        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto ci{inclusions[i]->volume_fraction()};
            const auto& S{dynamic_cast<eshelby_tensor_conductivity_base<_T, _Dim>*>(inclusions[i]->eshelby_tensor())->tensor()};
            const auto& Ci{make_conductivity_material(inclusions[i]->material())->conductivity_tensor()};
            TT[i] = tmech::inv(I + S*Cm_inv*(Ci - Cm));
            TT_sum += ci*TT[i];
        }
        TT_sum = tmech::eval(tmech::inv(TT_sum));

        for(size_type i{0}; i<inclusions.size(); ++i){
            A[i] = TT[i]*TT_sum;
        }
        A.back() = TT_sum;
    }
};

#endif // MEAN_FIELD_COMPOSITE_CONDUCTIVITY_KERNAL_BASE_BONES_H
