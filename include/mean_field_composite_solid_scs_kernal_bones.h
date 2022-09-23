#ifndef MEAN_FIELD_COMPOSITE_SOLID_SCS_KERNAL_BONES_H
#define MEAN_FIELD_COMPOSITE_SOLID_SCS_KERNAL_BONES_H

//self consistent
template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_scs_kernal :
        public mean_field_composite_solid_kernal_base<_T, _Dim, _Container>
{
public:
    using size_type  = std::size_t;
    using value_type = _T;

    mean_field_composite_scs_kernal() {}

    constexpr inline virtual void determine_strain_concentration_tensors(mean_field_composite_solid_base<_T, _Dim, _Container> & _composite) const override {
        auto& A{dynamic_cast<composite_material_solid_base<_T, _Dim, _Container>*>(&_composite)->strain_concentration_tensors()};
        const auto& inclusions{_composite.inclusions()};
        const auto* matrix_material{make_solid_material(_composite.matrix_material())};

        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};

        const auto& Cm{matrix_material->tangent_tensor()};
        const auto cm{_composite.matrix_volume_fraction()};

        std::vector<tmech::tensor<value_type, _Dim, 4>> Adil(inclusions.size());
        tmech::tensor<value_type, _Dim, 4> Adil_sum;
        tmech::tensor<value_type, _Dim, 4> Cm_inv{tmech::inv(Cm)};
        tmech::tensor<value_type, _Dim, 4> C;
        tmech::tensor<value_type, _Dim, 4> C_old;

        //start with mori-tanaka approximation
        Adil_sum = IIsym*cm;

        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto ci{inclusions[i]->volume_fraction()};
            const auto& S{dynamic_cast<eshelby_tensor_solid_base<_T, _Dim>*>(inclusions[i]->eshelby_tensor())->tensor()};
            const auto& Ci{make_solid_material(inclusions[i]->material())->tangent_tensor()};
            Adil[i] = tmech::inv(IIsym + tmech::dcontract(S, dcontract(Cm_inv, Ci - Cm)));
            Adil_sum += ci*Adil[i];
        }

        Adil_sum = tmech::eval(tmech::inv(Adil_sum));

        for(size_type i{0}; i<inclusions.size(); ++i){
            A[i] = tmech::dcontract(Adil[i], Adil_sum);
        }
        A.back() = tmech::dcontract(IIsym, Adil_sum);

        C = cm*tmech::dcontract(matrix_material->tangent_tensor(), A.back());

        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto ci{inclusions[i]->volume_fraction()};
            const auto& Ci{make_solid_material(inclusions[i]->material())->tangent_tensor()};
            C += ci*tmech::dcontract(Ci, A[i]);
        }

        size_type max_iter{100}, iter{0};
        while(true){
            const auto [K, mu]{extract_elasticity_parameter_K_mu(C)};
            const auto nu{nu_K_mu(K,mu)};
            const tmech::tensor<value_type, _Dim, 4> Cinv{tmech::inv(C)};
            for(size_type i{0}; i<inclusions.size(); ++i){
                const auto& Ci{make_solid_material(inclusions[i]->material())->tangent_tensor()};
                inclusions[i]->init();
                const auto& S{dynamic_cast<eshelby_tensor_solid_base<_T, _Dim>*>(inclusions[i]->eshelby_tensor())->tensor()};
                A[i] = 0.75*A[i] + 0.25*tmech::inv(IIsym + tmech::dcontract(S, tmech::dcontract(Cinv, (Ci - C))));
            }

            C_old = C;
            C = Cm;
            for(size_type i{0}; i<inclusions.size(); ++i){
                const auto ci{inclusions[i]->volume_fraction()};
                const auto& Ci{make_solid_material(inclusions[i]->material())->tangent_tensor()};
                C += ci*tmech::dcontract(Ci - Cm, A[i]);
            }

            const auto error{norm(C_old - C)};

            if(error <= 5e-5){
                break;
            }

            if(iter == max_iter){
                throw std::runtime_error("mean_field_composite_scs_kernal::determine_strain_concentration_tensors(): no convergenz");
            }
            ++iter;
        }

        tmech::tensor<value_type, _Dim, 4> Ascs_sum;
        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto ci{inclusions[i]->volume_fraction()};
            Ascs_sum += ci*A[i];
        }
        A.back() = (IIsym - Ascs_sum)/cm;
    }
};

#endif // MEAN_FIELD_COMPOSITE_SOLID_SCS_KERNAL_BONES_H
