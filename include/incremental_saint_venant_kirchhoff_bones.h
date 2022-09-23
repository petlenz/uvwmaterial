#ifndef INCREMENTAL_SAINT_VENANT_KIRCHHOFF_BONES_H
#define INCREMENTAL_SAINT_VENANT_KIRCHHOFF_BONES_H


template <typename _T, std::size_t _Dim, typename _Container>
class incremental_saint_venant_kirchhoff :
        public finite_strain_solid_material_base<_T, _Dim, _Container>,
        public incremental_solid_material<_T, _Dim>,
        public solid_material_base<_T, _Dim, _Container>

{
public:
    using value_type = _T;
    using size_type = std::size_t;

    incremental_saint_venant_kirchhoff():
        finite_strain_solid_material_base<_T, _Dim, _Container>(this)
    {}

    virtual ~incremental_saint_venant_kirchhoff(){}

    inline virtual void init(){
        if(!this->_is_init){
            this->_type = FINITE_STRAIN_FORMULATION::LAGRANGIAN;
            this->reinit();
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){
        const tmech::eye<value_type, _Dim, 2> I;
        this->_C = this->_parameter[0]*tmech::otimes(I,I) + this->_parameter[1]*(tmech::otimesu(I,I) + tmech::otimesl(I,I));
    }

    inline virtual void update(){
        const tmech::eye<value_type, _Dim, 2> I;
        //dE = E^{n+1} - E^{n}
        const tmech::tensor<value_type,_Dim,2> dE{0.5*(tmech::trans(this->_F)*this->_F - tmech::trans(this->_Fn)*this->_Fn)};
        this->_stress += tmech::dcontract(this->_C, dE);
    }

    inline virtual void update_stress(){
        const tmech::eye<value_type, _Dim, 2> I;
        const tmech::tensor<value_type,_Dim,2> dE{0.5*(tmech::trans(this->_F+this->_dstrain)*(this->_F+this->_dstrain) - tmech::trans(this->_Fn)*this->_Fn)};
        this->_stress += tmech::dcontract(this->_C, dE);
    }

    inline virtual void update_tangent(){
        //nothing to do
    }


    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override{
        return this;
    }
};

#endif // INCREMENTAL_SAINT_VENANT_KIRCHHOFF_BONES_H
