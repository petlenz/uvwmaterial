#ifndef THERMO_MECHANICAL_ELASTIC_BONES_H
#define THERMO_MECHANICAL_ELASTIC_BONES_H


//template<typename _T, std::size_t _Dim, typename _Container>
//class thermo_mechanical_elastic :
//        public thermo_mechanical_material_base<_T, _Dim, _Container>
//{
//public:
//    using value_type = _T;

//    thermo_mechanical_elastic() {}

//    virtual ~thermo_mechanical_elastic(){}

//    inline virtual void init(){
//        if(!this->_is_init){
//            const auto K{this->_parameter[0]};
//            const auto G{this->_parameter[1]};
//            const auto k{this->_parameter[3]};
//            constexpr _T fac{static_cast<_T>(_Dim)};
//            const auto I{tmech::eye<_T, _Dim, 2>()};
//            const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
//            const auto IIvol{tmech::otimes(I, I)/fac};
//            const auto IIdev{IIsym - IIvol};
//            this->_C = 3*K*IIvol + 2*G*IIdev;

//            this->_k = k*I;
//            this->_is_init = true;
//        }
//    }

//    inline virtual void reinit(){}

//    inline virtual void update(){
//        const tmech::eye<_T, _Dim, 2> I;
//        const auto alpha{this->_parameter[2]*I*this->_temperature};
//        this->_stress = tmech::dcontract(this->_C, this->_strain - alpha);
//        this->_q = this->_k*this->_dT;
//        //dq/du
//        //dq/dg
//        //dsigma/dT latent heat stress
//        //dsigma/deps
//    }

//    inline virtual void update_stress(){

//    }

//    inline virtual void update_tangent(){

//    }

//    inline virtual void update_heat_flux(){

//    }

//    inline virtual void update_conductivity(){

//    }

//    inline virtual thermo_mechanical_material_base<_T, _Dim, _Container>* base_material() override{
//        return this;
//    }
//};
#endif // THERMO_MECHANICAL_ELASTIC_BONES_H
