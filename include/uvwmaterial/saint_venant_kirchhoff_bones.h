#ifndef SAINT_VENANT_KIRCHHOFF_BONES_H
#define SAINT_VENANT_KIRCHHOFF_BONES_H


template <typename _T, std::size_t _Dim, typename _Container>
class saint_venant_kirchhoff :
        public solid_material_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>

{
public:
    using value_type = _T;
    using finite_strain_base = finite_strain_solid_material_base<_T, _Dim, _Container>;

    saint_venant_kirchhoff():
        solid_material_base<_T, _Dim, _Container>(),
        finite_strain_base(this)
    {}

    virtual ~saint_venant_kirchhoff(){}

    inline virtual void init(){
        if(!this->_is_init){
            this->_type = FINITE_STRAIN_FORMULATION::LAGRANGIAN;
            this->reinit();
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){
        const tmech::eye<value_type, _Dim, 2> I;
        this->_C = this->_parameter[0]*tmech::otimes(I,I) + this->_parameter[1]*(tmech::otimesu(I,I)+tmech::otimesl(I,I));
        //_C0 = this->_parameter[0]*tmech::otimes(I,I) + this->_parameter[1]*(tmech::otimesu(I,I)+tmech::otimesl(I,I));
    }

    inline virtual void update(){
        const tmech::eye<value_type, _Dim, 2> I;
        const auto E = 0.5*(tmech::trans(this->_F)*this->_F - I);
        this->_stress = tmech::dcontract(this->_C, E);
//        this->_C = tmech::dcontract(_C0,
//                                    (0.5*(tmech::basis_change<tmech::sequence<1,3,4,2>>(tmech::inner_product<tmech::sequence<2>,tmech::sequence<1>>(tmech::basis_change<tmech::sequence<2,1,3,4>>(tmech::otimesu(I, I)), this->_F))+tmech::inner_product<tmech::sequence<2>,tmech::sequence<1>>(tmech::trans(this->_F), tmech::otimesu(I, I)))))
//                + tmech::otimesu(I, this->_stress);
//        this->_stress = tmech::eval(this->_F*this->_stress);
    }

    inline virtual void update_stress(){
        const tmech::eye<value_type, _Dim, 2> I;
        const auto E = 0.5*(tmech::trans(this->_F)*this->_F - I);
        this->_stress = tmech::dcontract(this->_C, E);
    }

    inline virtual void update_tangent(){
        //nothing to do
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override{
        return this;
    }

private:
    tmech::tensor<value_type,_Dim,4> _C0;
};

#endif // SAINT_VENANT_KIRCHHOFF_BONES_H
