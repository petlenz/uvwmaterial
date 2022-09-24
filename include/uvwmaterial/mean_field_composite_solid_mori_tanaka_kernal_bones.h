/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef MEAN_FIELD_COMPOSITE_SOLID_MORI_TANAKA_KERNAL_BONES_H
#define MEAN_FIELD_COMPOSITE_SOLID_MORI_TANAKA_KERNAL_BONES_H

//mori tanaka
template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_mori_tanaka_kernal : public mean_field_composite_solid_kernal_base<_T, _Dim, _Container>
{
public:
    using size_type  = std::size_t;
    using value_type = _T;

    mean_field_composite_mori_tanaka_kernal() {}

    inline virtual void determine_strain_concentration_tensors(mean_field_composite_solid_base<_T, _Dim, _Container> & _composite) const override {
        auto& A{dynamic_cast<composite_material_solid_base<_T, _Dim, _Container>*>(&_composite)->strain_concentration_tensors()};
        const auto& inclusions{_composite.inclusions()};
        const auto* matrix_material{make_solid_material(_composite.matrix_material())};
        std::vector<tmech::tensor<_T, _Dim, 4>> TT(inclusions.size());

        tmech::tensor<value_type, _Dim, 4> TT_sum;
        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + tmech::otimesl(I,I))*0.5};

        //A_i = A_dil_i : inv(sum_r c_r A_dil_r)
        //r in [0, n+1]; i in [0, n+1]
        //A_dil_n+1 = cm*IIsym;

        TT_sum = IIsym*_composite.matrix_volume_fraction();

        const auto& Cm{matrix_material->tangent_tensor()};
        const tmech::tensor<value_type, _Dim, 4> Cm_inv{tmech::inv(Cm)};
        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto ci{inclusions[i]->volume_fraction()};
            const auto& S{dynamic_cast<eshelby_tensor_solid_base<_T, _Dim>*>(inclusions[i]->eshelby_tensor())->tensor()};
            const auto& Ci{make_solid_material(inclusions[i]->material())->tangent_tensor()};
            TT[i] = tmech::inv(IIsym +  tmech::dcontract(S, tmech::dcontract(Cm_inv, (Ci - Cm))));
            TT_sum += ci*TT[i];
        }
        TT_sum = tmech::eval(tmech::inv(TT_sum));

        for(size_type i{0}; i<inclusions.size(); ++i){
            A[i] = tmech::dcontract(TT[i], TT_sum);
        }
        A.back() = tmech::dcontract(IIsym, TT_sum);
    }
};

#endif // MEAN_FIELD_COMPOSITE_SOLID_MORI_TANAKA_KERNAL_BONES_H
