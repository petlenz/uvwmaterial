#ifndef FINITE_STRAIN_SOLID_MATERIAL_BASE_BONES_H
#define FINITE_STRAIN_SOLID_MATERIAL_BASE_BONES_H

enum FINITE_STRAIN_FORMULATION{LAGRANGIAN, MIXED, SPATIAL};


template <typename _T, std::size_t _Dim, typename _Container>
class finite_strain_solid_material_base //:
        //public solid_material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using tensor2 = tmech::tensor<value_type, _Dim, 2>;

    finite_strain_solid_material_base(material_base<_T, _Dim, _Container>* __solid_material):
        _F(tmech::eye<_T,_Dim,2>()),
        _Fn(tmech::eye<_T,_Dim,2>()),
        _type(FINITE_STRAIN_FORMULATION::MIXED),
        _solid_material(dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(__solid_material)),
        _thermal_material(dynamic_cast<thermo_mechanical_material_base<_T, _Dim, _Container>*>(__solid_material))
    {}

    virtual ~finite_strain_solid_material_base(){}

    constexpr inline auto const& deformation_tensor()const{
        return _F;
    }

    constexpr inline auto& deformation_tensor(){
        return _F;
    }

    constexpr inline auto const& deformation_tensor_n()const{
        return _Fn;
    }

    constexpr inline auto& deformation_tensor_n(){
        return _Fn;
    }

    constexpr inline auto formulation_type()const{
        return _type;
    }

    constexpr inline tmech::tensor<value_type,_Dim,4> mixed_tangent_tensor()const{
        const tmech::eye<value_type,_Dim,2> I;
        switch (_type) {
        case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
            //return  tmech::dcontract(_solid_material->tangent_tensor(),
            //                         (0.5*(tmech::basis_change<tmech::sequence<1,3,4,2>>(tmech::inner_product<tmech::sequence<2>,tmech::sequence<1>>(tmech::basis_change<tmech::sequence<2,1,3,4>>(tmech::otimesu(I, I)), _F))+tmech::inner_product<tmech::sequence<2>,tmech::sequence<1>>(tmech::trans(_F), tmech::otimesu(I, I)))))
            //        + tmech::otimesu(I, _solid_material->stress_tensor());
            return tmech::dcontract(tmech::dcontract(tmech::otimesu(_F, I), _solid_material->tangent_tensor()), tmech::otimesu(tmech::trans(_F), I))
                    + tmech::otimesu(I, _solid_material->stress_tensor());
            break;
        default:
            return _solid_material->tangent_tensor();
            break;
        }
    }

    constexpr inline tmech::tensor<value_type,_Dim,2> mixed_stress_tensor()const{
        tmech::tensor<value_type,_Dim,2> temp;
        switch (_type) {
        case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
            temp = _F*_solid_material->stress_tensor();
            break;
        default:
            temp = _solid_material->stress_tensor();
            break;
        }
        return temp;
    }


    constexpr inline tmech::tensor<value_type,_Dim,4> lagrangian_tangent_tensor()const{
        const tmech::eye<value_type,_Dim,2> I;
        switch (_type) {
        case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
            //return  tmech::dcontract(_solid_material->tangent_tensor(),
            //                         (0.5*(tmech::basis_change<tmech::sequence<1,3,4,2>>(tmech::inner_product<tmech::sequence<2>,tmech::sequence<1>>(tmech::basis_change<tmech::sequence<2,1,3,4>>(tmech::otimesu(I, I)), _F))+tmech::inner_product<tmech::sequence<2>,tmech::sequence<1>>(tmech::trans(_F), tmech::otimesu(I, I)))))
            //        + tmech::otimesu(I, _solid_material->stress_tensor());
            return _solid_material->tangent_tensor();
            break;
        default:
            return tmech::dcontract(tmech::dcontract(tmech::invf(tmech::otimesu(_F, I)), _solid_material->tangent_tensor() - tmech::otimesu(I, _solid_material->stress_tensor())), tmech::invf(tmech::otimesu(tmech::trans(_F), I)));
            break;
        }
    }

    constexpr inline tmech::tensor<value_type,_Dim,2> lagrangian_stress_tensor()const{
        switch (_type) {
        case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
            return _solid_material->stress_tensor();
            break;
        default:
            return tmech::inv(_F)*_solid_material->stress_tensor();
            break;
        }
    }

    constexpr inline tmech::tensor<value_type,_Dim,4> spatial_tangent_tensor()const{
        const tmech::eye<value_type,_Dim,2> I;
        switch (_type) {
        case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
            return (tmech::dcontract(tmech::otimesu(_F,_F),tmech::dcontract(_solid_material->tangent_tensor(), tmech::otimesu(tmech::trans(_F),tmech::trans(_F)))))/tmech::det(_F);
        case FINITE_STRAIN_FORMULATION::MIXED:
            return tmech::dcontract(tmech::dcontract(tmech::invf(tmech::otimesu(_F, I)), _solid_material->tangent_tensor() - tmech::otimesu(I, _solid_material->stress_tensor())), tmech::invf(tmech::otimesu(tmech::trans(_F), I)));
        default:
            return _solid_material->tangent_tensor();
        }
    }


    constexpr inline tmech::tensor<value_type,_Dim,2> spatial_stress_tensor()const{
        if(_solid_material){
            switch (_type) {
            case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
                return _F*_solid_material->stress_tensor()*tmech::trans(_F)/tmech::det(_F);
            case FINITE_STRAIN_FORMULATION::MIXED:
                return _solid_material->stress_tensor()*tmech::trans(_F)/tmech::det(_F);
            default:
                return _solid_material->stress_tensor();
            }
        }else{
            switch (_type) {
            case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
                return _F*_thermal_material->stress_tensor()*tmech::trans(_F)/tmech::det(_F);
            case FINITE_STRAIN_FORMULATION::MIXED:
                return _thermal_material->stress_tensor()*tmech::trans(_F)/tmech::det(_F);
            default:
                return _thermal_material->stress_tensor();
            }
        }
    }

protected:

    tensor2 _F;
    tensor2 _Fn;
    FINITE_STRAIN_FORMULATION _type;
    solid_material_base<_T, _Dim, _Container>* _solid_material;
    thermo_mechanical_material_base<_T, _Dim, _Container>* _thermal_material;

};

#endif // FINITE_STRAIN_SOLID_MATERIAL_BASE_BONES_H
