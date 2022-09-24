/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef MEAN_FIELD_COMPOSITE_SOLID_FINITE_STRAIN_DILUTE_KERNAL_BONES_H
#define MEAN_FIELD_COMPOSITE_SOLID_FINITE_STRAIN_DILUTE_KERNAL_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_solid_finite_strain_dilute_kernal :
        public mean_field_composite_solid_kernal_base<_T, _Dim, _Container>
{
public:
    using size_type  = std::size_t;
    using value_type = _T;

    mean_field_composite_solid_finite_strain_dilute_kernal() {}

    inline virtual void determine_strain_concentration_tensors(mean_field_composite_solid_base<_T, _Dim, _Container> & _composite) const override {
        auto& A{dynamic_cast<composite_material_solid_base<_T, _Dim, _Container>*>(&_composite)->strain_concentration_tensors()};
        const auto& inclusions{_composite.inclusions()};
        const auto* matrix_material{make_solid_material(_composite.matrix_material())};

        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto II{tmech::otimesu(I,I)};
        const auto& Cm{matrix_material->tangent_tensor()};
        A.back() = II;

        const tmech::tensor<value_type, _Dim, 4> Cm_inv{tmech::inv(Cm)};
        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto& Ci{make_solid_material(inclusions[i]->material())->tangent_tensor()};
            const auto& S{dynamic_cast<eshelby_tensor_solid_base<_T, _Dim>*>(inclusions[i]->eshelby_tensor())->tensor()};
            A[i] = tmech::inv(II + tmech::dcontract(S, tmech::dcontract(Cm_inv, Ci - Cm)));
        }
    }
};


template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_solid_finite_strain_mori_tanaka_kernal :
        public mean_field_composite_solid_kernal_base<_T, _Dim, _Container>
{
public:
    using size_type  = std::size_t;
    using value_type = _T;

    mean_field_composite_solid_finite_strain_mori_tanaka_kernal() {}

    constexpr inline virtual void determine_strain_concentration_tensors(mean_field_composite_solid_base<_T, _Dim, _Container> & _composite) const override {
        auto& A{dynamic_cast<composite_material_solid_base<_T, _Dim, _Container>*>(&_composite)->strain_concentration_tensors()};
        const auto& inclusions{_composite.inclusions()};
        const auto* matrix_material{dynamic_cast<const finite_strain_solid_material_base<_T,_Dim,_Container>*>(_composite.matrix_material())};
        //const auto* matrix_material{dynamic_cast<const solid_material_base<_T,_Dim,_Container>*>(_composite.matrix_material())};
        std::vector<tmech::tensor<_T, _Dim, 4>> TT(inclusions.size());

        tmech::tensor<value_type, _Dim, 4> TT_sum;
        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto II{tmech::otimesu(I,I)};

        //A_i = A_dil_i : inv(sum_r c_r A_dil_r)
        //r in [0, n+1]; i in [0, n+1]
        //A_dil_n+1 = cm*IIsym;

        TT_sum = II*_composite.matrix_volume_fraction();
        //std::cout<<"matrix "<<_composite.matrix_volume_fraction()<<std::endl;

        //tmech::tensor<_T, _Dim, 4> Cm{matrix_material->tangent_tensor()};
        tmech::tensor<_T, _Dim, 4> Cm{matrix_material->mixed_tangent_tensor()};
        tmech::tensor<value_type, _Dim, 4> Cm_inv;

        auto is_minor_sym = [](auto const& A){
            constexpr auto Dim{A.dimension()};
            for(std::size_t i{0}; i<Dim; ++i){
                for(std::size_t j{0}; j<Dim; ++j){
                    for(std::size_t k{0}; k<Dim; ++k){
                        for(std::size_t l{0}; l<Dim; ++l){
                            if(std::abs(A(i,j,k,l) - A(j,i,k,l)) > std::numeric_limits<value_type>::epsilon()*std::max((_T)1.0, std::max(A(i,j,k,l), A(j,i,k,l)))*10.0 ||
                                    std::abs(A(i,j,k,l) - A(i,j,l,k)) > std::numeric_limits<value_type>::epsilon()*std::max((_T)1.0, std::max(A(i,j,k,l), A(i,j,l,k)))*10.0 ){
                                //std::cout<<std::abs(A(i,j,k,l) - A(j,i,k,l))<<" "<<std::abs(A(i,j,k,l) - A(i,j,l,k))<<" "<<std::numeric_limits<value_type>::epsilon()*std::max(A(i,j,k,l), A(j,i,k,l))<<std::endl;
                                return false;
                            }
                        }
                    }
                }
            }
            return true;
        };
        if(is_minor_sym(Cm)){
            Cm_inv = tmech::inv(Cm);
        }else{
            Cm_inv = tmech::invf(Cm);
        }
        //std::cout<<"matrix \n"<<Cm<<std::endl;
        //std::cout<<"Mixed\n";
        //std::cout<<tmech::dcontract(Cm, matrix_material->deformation_tensor())<<std::endl;
        //std::cout<<matrix_material->mixed_stress_tensor()<<std::endl;
        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto ci{inclusions[i]->volume_fraction()};
            //std::cout<<"inclusion "<<ci<<std::endl;

            const auto& S{dynamic_cast<eshelby_tensor_solid_base<_T, _Dim>*>(inclusions[i]->eshelby_tensor())->tensor()};
            tmech::tensor<_T, _Dim, 4> Ci{dynamic_cast<const finite_strain_solid_material_base<_T,_Dim,_Container>*>(inclusions[i]->material())->mixed_tangent_tensor()};
            //std::cout<<"inclusion \n"<<Ci<<std::endl;
            //const auto& Ci{dynamic_cast<const solid_material_base<_T,_Dim,_Container>*>(inclusions[i]->material())->tangent_tensor()};
            TT[i] = tmech::invf(II +  tmech::dcontract(S, tmech::dcontract(Cm_inv, (Ci - Cm))));
            //std::cout<<S<<std::endl;
            TT_sum += ci*TT[i];
        }
        TT_sum = tmech::eval(tmech::invf(TT_sum));

        for(size_type i{0}; i<inclusions.size(); ++i){
            A[i] = tmech::dcontract(TT[i], TT_sum);
        }
        A.back() = tmech::dcontract(II, TT_sum);

        //std::cout<<0.7*A[0] + 0.3*A[1]<<std::endl;
        //std::cout<<tmech::otimesu(I,I)<<std::endl;
    }
};


#endif // MEAN_FIELD_COMPOSITE_SOLID_FINITE_STRAIN_DILUTE_KERNAL_BONES_H
