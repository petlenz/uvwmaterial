/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef COMPOSITE_MATERIAL_RULE_OF_MIXTURE_KERNAL_BASE_MEAT_H
#define COMPOSITE_MATERIAL_RULE_OF_MIXTURE_KERNAL_BASE_MEAT_H

template <typename _T, std::size_t _Dim, typename _Container>
rule_of_mixture_voigt<_T, _Dim, _Container>::rule_of_mixture_voigt(){}

template <typename _T, std::size_t _Dim, typename _Container>
rule_of_mixture_voigt<_T, _Dim, _Container>::~rule_of_mixture_voigt(){}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline void rule_of_mixture_voigt<_T, _Dim, _Container>::determine_strain_concentration_tensors(
        std::vector<small_strain_material_base<_T, _Dim, _Container>*> const& /*_materials*/,
        std::vector<_T> const& /*_volume_fractions*/,
        std::vector<tmech::tensor<_T, _Dim, 4>> & _A) const {

    const auto I{tmech::eye<value_type, _Dim, 2>()};
    const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
    for(size_type i{0}; i<_A.size(); ++i){
        _A[i] = IIsym;
    }
}



template <typename _T, std::size_t _Dim, typename _Container>
rule_of_mixture_reuss<_T, _Dim, _Container>::rule_of_mixture_reuss(){}

template <typename _T, std::size_t _Dim, typename _Container>
rule_of_mixture_reuss<_T, _Dim, _Container>::~rule_of_mixture_reuss(){}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline void rule_of_mixture_reuss<_T, _Dim, _Container>::determine_strain_concentration_tensors(std::vector<small_strain_material_base<_T, _Dim, _Container>*> const& _materials,
                                                                     std::vector<_T> const& _volume_fractions,
                                                                     std::vector<tmech::tensor<_T, _Dim, 4>> & _A) const {
    tmech::tensor<value_type, _Dim, 4> C;
    for(size_type i{0}; i<_materials.size(); ++i){
        _A[i] = tmech::inv(_materials[i]->tangent_tensor());
        C += _volume_fractions[i]*_A[i];
    }

    C = tmech::eval(tmech::inv(C));

    for(size_type i{0}; i<_materials.size(); ++i){
        _A[i] =  tmech::eval(tmech::dcontract(_A[i], C));
    }
}


template <typename _T, std::size_t _Dim, typename _Container>
rule_of_mixture_vrh<_T, _Dim, _Container>::rule_of_mixture_vrh(){}

template <typename _T, std::size_t _Dim, typename _Container>
rule_of_mixture_vrh<_T, _Dim, _Container>::~rule_of_mixture_vrh(){}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline void rule_of_mixture_vrh<_T, _Dim, _Container>::determine_strain_concentration_tensors(std::vector<small_strain_material_base<_T, _Dim, _Container>*> const& _materials,
                                                                                                        std::vector<_T> const& _volume_fractions,
                                                                                                        std::vector<tmech::tensor<_T, _Dim, 4>> & _A) const {
    const auto I{tmech::eye<value_type, _Dim, 2>()};
    const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
    tmech::tensor<value_type, _Dim, 4> C_R;
    C_R.fill(0);

    //A_i^V = IIsym
    //A_i^R = inv(C_i):inv(sum_i c_i inv(C_i))
    //A_i = 0.5*(A_i^V + A_i^R)
    //C = sum c_i (C_i + inv(C_i))
    for(size_type i{0}; i<_materials.size(); ++i){
        _A[i] = tmech::inv(_materials[i]->tangent_tensor());
        C_R += _volume_fractions[i]*_A[i];
    }

    C_R = tmech::eval(tmech::inv(C_R));

    //setup strainconcentration tensors
    for(size_type i{0}; i<_materials.size(); ++i){
        _A[i] = (IIsym + tmech::dcontract(_A[i], C_R))*0.5;
    }
}
#endif // COMPOSITE_MATERIAL_RULE_OF_MIXTURE_KERNAL_BASE_MEAT_H
