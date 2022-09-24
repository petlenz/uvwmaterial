/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef POLYMER_MATRIX_CONDUCTIVITY_MATERIAL_BONES_H
#define POLYMER_MATRIX_CONDUCTIVITY_MATERIAL_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class polymer_matrix_conductivity_material :
        public conductivity_material_base<_T, _Dim, _Container>,
        public history_material_base<_T>,
        public time_dependent_material_base<_T>
{
    enum CompositeType{MEAN_FIELD, RULE_OF_MIXTURE};
public:
    polymer_matrix_conductivity_material():
        conductivity_material_base<_T, _Dim, _Container>(),
        history_material_base<_T>(1)
    {}

    inline virtual void init(){
        if(!this->_is_init){
            if(dynamic_cast<mean_field_composite_conductivity_base<_T, _Dim, _Container>*>(_composite)){
                _composite_type = CompositeType::MEAN_FIELD;
            }/*else if(dynamic_cast<>()){
                _composite_type = CompositeType::RULE_OF_MIXTURE;
            }*/else{
                throw std::runtime_error("polymer_matrix_conductivity_material::init() no matching composite material");
            }
            dynamic_cast<conductivity_material_base<_T, _Dim, _Container>*>(_composite)->init();
            this->_is_init = true;
        }
    }

    inline virtual void reinit()override{
        dynamic_cast<conductivity_material_base<_T, _Dim, _Container>*>(_composite)->reinit();
    }

    inline virtual void update(){
        //update curing
        _curing_function->set_temperature(this->_scalar);
        _curing_function->set_history(this->_history[0]);
        _curing_function->set_time_increment(this->_delta_time);
        _curing_function->update();
        this->_history[0] = _curing_function->get_degree_of_cure();

        //update volume fractions
        const auto n{this->_parameter[0]};
        const auto cs{this->_history[0]};
        const auto cc{(1-n)*(1-cs)};
        const auto cr{n*(1-cs)};


        _composite->volume_fraction(0) = 1.0 - cr; //inclusion
        _composite->volume_fraction(1) = cr; //matrix
        //get inclusion
        auto sc_composite = dynamic_cast<composite_material_conductivity_base<_T, _Dim, _Container>*>(_composite->material(0));
        sc_composite->volume_fraction(0) = cs/(cc + cr); //inclusion
        sc_composite->volume_fraction(1) = cs + cr; //matrix
        this->reinit();
        auto mat = dynamic_cast<conductivity_material_base<_T, _Dim, _Container>*>(_composite);
        mat->update();
        this->_flux = mat->flux_tensor();
        this->_conductivity = mat->conductivity_tensor();
        this->_history[1] = this->_conductivity(0,0);
    }

    inline virtual void update_flux(){

    }

    inline virtual void update_conductivity(){

    }

    inline virtual conductivity_material_base<_T, _Dim, _Container>* base_material() {
        return this;
    }

    constexpr inline auto set_curing_function(kinetic_model_thermoset_base<_T>* __curing_function){
        _curing_function = __curing_function;
    }

    constexpr inline auto set_composite_material(composite_material_conductivity_base<_T, _Dim, _Container>* __composite){
        _composite = __composite;
    }

private:
    constexpr inline auto update_volume_fraction(){
        if(_composite_type == CompositeType::MEAN_FIELD){
            //_composite->volume_fraction(0) = 1.0 - ;
            auto sc_composite = dynamic_cast<composite_material_conductivity_base<_T, _Dim, _Container>*>(_composite->material(0));
//            auto solid = sc_composite->material(0);
//            auto curing = sc_composite->material(1);
//            auto resin = _composite->material(1);
        }
    }

protected:
    kinetic_model_thermoset_base<_T>* _curing_function;
    composite_material_conductivity_base<_T, _Dim, _Container>* _composite;
    CompositeType _composite_type;
};

#endif // POLYMER_MATRIX_CONDUCTIVITY_MATERIAL_BONES_H
