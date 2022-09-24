/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef HENCKY_STRAINS_LINEAR_ELASTIC_BONES_H
#define HENCKY_STRAINS_LINEAR_ELASTIC_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class hencky_strains_linear_elastic :
        public solid_material_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>

{
public:
    using value_type = _T;
    using finite_strain_base = finite_strain_solid_material_base<_T, _Dim, _Container>;

    hencky_strains_linear_elastic():
        solid_material_base<_T, _Dim, _Container>(),
        finite_strain_base(this)
    {}

    virtual ~hencky_strains_linear_elastic(){}

    inline virtual void init(){
        if(!this->_is_init){
            this->_type = FINITE_STRAIN_FORMULATION::LAGRANGIAN;
            this->reinit();
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){
        constexpr value_type fac{static_cast<value_type>(_Dim)};
        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
        const auto IIvol{tmech::otimes(I, I)/fac};
        const auto IIdev{IIsym - IIvol};
        const auto _K{this->_parameter[0]};
        const auto _G{this->_parameter[1]};
        _C0 = 3*_K*IIvol + 2*_G*IIdev;
        this->_C = _C0;
    }

    inline virtual void update(){
        auto logE = tmech::log(tmech::trans(this->_F)*this->_F);
        const auto E = 0.5*logE;
        const auto P = logE.derivative();
        const auto L = 2*logE.second_derivative();
        const auto T = tmech::dcontract(_C0, E);
        this->_stress = tmech::dcontract(P, T);
        this->_C = tmech::dcontract(P,tmech::dcontract(_C0,P)) + tmech::dcontract(T,L);
    }

    inline virtual void update_stress(){
        auto logE = tmech::log(tmech::trans(this->_F)*this->_F);
        const auto E = 0.5*logE;
        const auto P = logE.derivative();
        const auto T = tmech::dcontract(_C0, E);
        this->_stress = tmech::dcontract(P, T);
    }

    inline virtual void update_tangent(){
        update();
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override{
        return this;
    }

private:
    tmech::tensor<value_type,_Dim,4> _C0;
};



template <typename _T, std::size_t _Dim, typename _Container>
class hencky_strains_linear_thermo_elastic :
        public solid_material_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>,
        public temperature_dependent_material_base<_T>

{
public:
    using value_type = _T;
    using finite_strain_base = finite_strain_solid_material_base<_T, _Dim, _Container>;

    hencky_strains_linear_thermo_elastic():
        solid_material_base<_T, _Dim, _Container>(),
        finite_strain_base(this),
        temperature_dependent_material_base<_T>()
    {}

    virtual ~hencky_strains_linear_thermo_elastic(){}

    inline virtual void init(){
        if(!this->_is_init){
            this->_type = FINITE_STRAIN_FORMULATION::LAGRANGIAN;
            this->reinit();
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){
        constexpr value_type fac{static_cast<value_type>(_Dim)};
        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
        const auto IIvol{tmech::otimes(I, I)/fac};
        const auto IIdev{IIsym - IIvol};
        const auto _K{this->_parameter[0]};
        const auto _G{this->_parameter[1]};
        _C0 = 3*_K*IIvol + 2*_G*IIdev;
        this->_C = _C0;
    }

    inline virtual void update(){
        const tmech::eye<_T,_Dim,2> I;
        auto logE = tmech::log(tmech::trans(this->_F)*this->_F);
        const auto E = 0.5*logE;
        const auto P = logE.derivative();
        const auto L = 2*logE.second_derivative();
        const auto Ea = this->_parameter[2]*this->_temperature*I;
        const auto T = tmech::dcontract(_C0, E - Ea);
        this->_stress = tmech::dcontract(P, T);
        this->_C = tmech::dcontract(P,tmech::dcontract(_C0,P)) + tmech::dcontract(T,L);
    }

    inline virtual void update_stress(){
        const tmech::eye<_T,_Dim,2> I;
        auto logE = tmech::log(tmech::trans(this->_F)*this->_F);
        const auto E = 0.5*logE;
        const auto P = logE.derivative();
        const auto Ea = this->_parameter[2]*this->_temperature*I;
        const auto T = tmech::dcontract(_C0, E - Ea);
        this->_stress = tmech::dcontract(P, T);
    }

    inline virtual void update_tangent(){
        update();
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override{
        return this;
    }

private:
    tmech::tensor<value_type,_Dim,4> _C0;
};



template <typename _T, std::size_t _Dim, typename _Container>
class hencky_strains_linear_thermo_chemo_elastic :
        public solid_material_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>,
        public temperature_dependent_material_base<_T>,
        public degree_of_cure_dependent_material_base<_T>
{
public:
    using value_type = _T;
    using finite_strain_base = finite_strain_solid_material_base<_T, _Dim, _Container>;

    hencky_strains_linear_thermo_chemo_elastic():
        solid_material_base<_T, _Dim, _Container>(),
        finite_strain_base(this),
        temperature_dependent_material_base<_T>(),
        degree_of_cure_dependent_material_base<_T>()
    {}

    virtual ~hencky_strains_linear_thermo_chemo_elastic(){}

    inline virtual void init(){
        if(!this->_is_init){
            this->_type = FINITE_STRAIN_FORMULATION::LAGRANGIAN;
            this->reinit();
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){
        constexpr value_type fac{static_cast<value_type>(_Dim)};
        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
        const auto IIvol{tmech::otimes(I, I)/fac};
        const auto IIdev{IIsym - IIvol};
        const auto _K{this->_parameter[0]};
        const auto _G{this->_parameter[1]};
        _C0 = 3*_K*IIvol + 2*_G*IIdev;
        this->_C = _C0;
    }

    inline virtual void update(){
        const tmech::eye<_T,_Dim,2> I;
        auto logE = tmech::log(tmech::trans(this->_F)*this->_F);
        const auto E = 0.5*logE;
        const auto P = logE.derivative();
        const auto L = 2*logE.second_derivative();
        const auto Ea = this->_parameter[2]*this->_temperature*I;
        const auto Eb = this->_parameter[3]*this->_curing*I;
        const auto T = tmech::dcontract(_C0, E - Ea - Eb);
        this->_stress = tmech::dcontract(P, T);
        this->_C = tmech::dcontract(P,tmech::dcontract(_C0,P)) + tmech::dcontract(T,L);
    }

    inline virtual void update_stress(){
        const tmech::eye<_T,_Dim,2> I;
        auto logE = tmech::log(tmech::trans(this->_F)*this->_F);
        const auto E = 0.5*logE;
        const auto P = logE.derivative();
        const auto Ea = this->_parameter[2]*this->_temperature*I;
        const auto Eb = this->_parameter[3]*this->_dcuring*I;
        const auto T = tmech::dcontract(_C0, E - Ea - Eb);
        this->_stress = tmech::dcontract(P, T);
    }

    inline virtual void update_tangent(){
        update();
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override{
        return this;
    }

private:
    tmech::tensor<value_type,_Dim,4> _C0;
};

#endif // HENCKY_STRAINS_LINEAR_ELASTIC_BONES_H
