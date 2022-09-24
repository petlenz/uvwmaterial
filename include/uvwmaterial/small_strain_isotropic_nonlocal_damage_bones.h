/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SMALL_STRAIN_ISOTROPIC_NONLOCAL_DAMAGE_BONES_H
#define SMALL_STRAIN_ISOTROPIC_NONLOCAL_DAMAGE_BONES_H

template <typename T, std::size_t Dim, typename Container>
class small_strain_isotropic_nonlocal_damage :
        public history_material_base<T>,
        public small_strain_material_base<T, Dim, Container>,
        public nonlocal_material_base<T, Dim>
{
public:
    using value_type = T;
    using size_type = std::size_t;
    using tensor4 = tmech::tensor<value_type, Dim, 4>;

    small_strain_isotropic_nonlocal_damage():
        history_material_base<T>(3),
        small_strain_material_base<T, Dim, Container>(),
        nonlocal_material_base<T, Dim>(),
        _material(nullptr),
        _damage_func(nullptr),
        _yield_function(nullptr)
    {}

    small_strain_isotropic_nonlocal_damage(small_strain_material_base<T, Dim, Container> * material, propagation_law_base<T> * damage_func, yield_function_isotropic_damage_strain_based<T, Dim, Container> * yield_function):
        history_material_base<T>(3),
        small_strain_material_base<T, Dim, Container>(),
        nonlocal_material_base<T, Dim>(),
        _material(material),
        _damage_func(damage_func),
        _yield_function(yield_function)
    {}

    virtual ~small_strain_isotropic_nonlocal_damage(){}

    inline virtual void init()override{
        if(!this->_is_init){
            _material->init();
            _ptr_nonlocal_var = &this->_nonlocal_variables.find(this->_nonlocal_domains[0])->second[0];
            //this->set_nonlocal_domain(this->_nonlocal_domains[0], 1);
            this->_C = _material->tangent_tensor();
            this->_is_init = true;
        }
    }

    inline virtual void reinit()override{
        _material->reinit();
        this->_C = _material->tangent_tensor();
    }

    inline virtual void update()override{
        //update base material
        //compute trail stress
        _material->strain_tensor() = this->_strain;
        _material->update();

        if(_average_type == "equivalent_strain"){
            update_equivalent_strain();
        }else if(_average_type == "damage_variable"){
            update_damage_variable();
        }else {
            throw std::runtime_error("small_strain_isotropic_nonlocal_damage::update(): no matching average type");
        }
    }

    inline virtual void update_stress()override{
        throw std::runtime_error("small_strain_isotropic_nonlocal_d_damage::update_stress()");
        if(_average_type == "equivalent_strain"){
            update_stress_equivalent_strain();
        }else if(_average_type == "damage_variable"){
            update_stress_damage_variable();
        }else {
            throw std::runtime_error("small_strain_isotropic_nonlocal_damage::update(): no matching average type");
        }
    }

    inline virtual void update_tangent()override{
        throw std::runtime_error("small_strain_isotropic_nonlocal_d_damage::update_tangent()");
        if(_average_type == "equivalent_strain"){
            update_tangent_equivalent_strain();
        }else if(_average_type == "damage_variable"){
            update_tangent_damage_variable();
        }else {
            throw std::runtime_error("small_strain_isotropic_nonlocal_damage::update(): no matching average type");
        }
    }

    constexpr inline auto set_base_material(small_strain_material_base<value_type, Dim, Container> * __material){
        _material = __material;
    }

    constexpr inline auto set_propagation_law(propagation_law_base<value_type> * __propagation_law){
        _damage_func = __propagation_law;
    }

    constexpr inline auto set_yield_function(yield_function_isotropic_damage_strain_based<value_type, Dim, Container> * __yield_function){
        _yield_function = __yield_function;
    }

    inline auto& average_type(){
        return _average_type;
    }

    inline auto const& average_type()const{
        return _average_type;
    }

    virtual inline solid_material_base<T, Dim, Container>* base_material(){
        return _material->base_material();
    }

    inline void update_nonlocal_variables() override {
        if(_average_type == "equivalent_strain"){
            update_nonlocal_variable_equivalent_strain();
            return ;
        }
        if(_average_type == "damage_variable"){
            update_nonlocal_variable_damage_variable();
            return ;
        }
        throw std::runtime_error("small_strain_isotropic_nonlocal_damage::update(): no matching average type");
    }

protected:
    std::string _average_type;
    small_strain_material_base<value_type, Dim, Container> * _material;
    propagation_law_base<value_type> * _damage_func;
    yield_function_isotropic_damage_strain_based<value_type, Dim, Container> * _yield_function;

private:
    inline auto update_equivalent_strain(){
        auto& sig{this->_stress};
        auto& C{this->_C};
        auto& history{this->_history};
        const auto& nonlocal_domain{this->_nonlocal_domains[0]};
        const auto nonlocal_val{*_ptr_nonlocal_var};
        auto& source{this->_source.find(nonlocal_domain)->second[0]};
        auto& receiver{this->_receiver.find(nonlocal_domain)->second[0]};
        auto& marker{this->_element_marker.find(nonlocal_domain)->second};

        if(nonlocal_val > std::max(history[1], _yield_function->critical_value())){
            _yield_function->update_equivalent_scalar();
            _yield_function->update_strain_dependent_parts();
            const auto D{_damage_func->value(nonlocal_val)};
            const auto dD{_damage_func->derivative(nonlocal_val)};
            C = (1.0-D)*_material->tangent_tensor();
            sig = (1.0-D)*_material->stress_tensor();
            receiver = -dD*_material->stress_tensor();
            source =  _yield_function->derivative();
            history[0] = D;
            history[1] = nonlocal_val;
            marker = true;
        }else{
            //No damage
            const auto D{history[0]};
            C = (1-D)*_material->tangent_tensor();
            sig = (1-D)*_material->stress_tensor();
            receiver.fill(0);
            source.fill(0);
            marker = false;
            history[1] = std::max(_yield_function->critical_value(), history[1]);
        }
        if(history[0] > 1.0){
            throw std::runtime_error("Nonlocal equivalent strain damage D>1");
        }
    }

    inline auto update_damage_variable(){
        auto& sig{this->_stress};
        auto& C{this->_C};
        auto& history{this->_history};
        const auto& nonlocal_domain{this->_nonlocal_domains[0]};
        const auto nonlocal_val{*_ptr_nonlocal_var};
        auto& source{this->source_tensor(nonlocal_domain)[0]};
        auto& receiver{this->receiver_tensor(nonlocal_domain)[0]};
        auto& marker{this->_element_marker.find(nonlocal_domain)->second};

        //update equivalent scalar
        _yield_function->update_equivalent_scalar();

        //solve for new damage multiplier
        //const auto his_new{_yield_function->solve(history[1])};
        //std::cout<<"His new "<<his_new<<" his old "<<history[1]<<std::endl;
        if(/*_yield_function->yielding(history[1])*/ _yield_function->state() >= std::max(_yield_function->critical_value(), history[1])){
            _yield_function->update_strain_dependent_parts();
            const auto D{nonlocal_val};//std::max(nonlocal_val, history[0])};
            const auto dD{_damage_func->derivative(_yield_function->state())};
            C = (1.0-D)*_material->tangent_tensor();
            sig = (1.0-D)*_material->stress_tensor();
            receiver = -_material->stress_tensor();
            source = dD*_yield_function->derivative();
            history[0] = _damage_func->value(_yield_function->state());
            history[2] = D;
            history[1] = _yield_function->state();
            //history[2] = _damage_func->value(_yield_function->state());
            marker = true;
        }else{
            //No damage
            const auto D{history[2]};
            C = (1-D)*_material->tangent_tensor();
            sig = (1-D)*_material->stress_tensor();
            receiver.fill(0);
            source.fill(0);
            marker = false;
            history[1] = std::max(_yield_function->critical_value(), history[1]);
        }

        if(history[0] > 1.0){
            throw std::runtime_error("Nonlocal damage variable damage D>1");
        }
    }

    inline auto update_tangent_equivalent_strain(){

    }

    inline auto update_tangent_damage_variable(){

    }

    inline auto update_stress_equivalent_strain(){

    }

    inline auto update_stress_damage_variable(){

    }


    constexpr inline auto update_nonlocal_variable_equivalent_strain(){
        //const auto& nonlocal_domain{this->_nonlocal_domains[0]};
        //auto& nonlocal_val{this->_nonlocal_variables.find(nonlocal_domain)->second[0]};
        _material->strain_tensor() = this->_strain;
        _material->update();
        _yield_function->update_equivalent_scalar();
        *_ptr_nonlocal_var = _yield_function->state();
    }

    constexpr inline auto update_nonlocal_variable_damage_variable(){
        //const auto& nonlocal_domain{this->_nonlocal_domains[0]};
        //auto& nonlocal_val{this->_nonlocal_variables.find(nonlocal_domain)->second[0]};
        _material->strain_tensor() = this->_strain;
        _material->update();
        _yield_function->update_equivalent_scalar();
        if(_yield_function->state() >= std::max(_yield_function->critical_value(), this->_history[1])){
            *_ptr_nonlocal_var = _damage_func->value(_yield_function->state());
        }else {
            *_ptr_nonlocal_var = this->_history[0];
        }
    }

    value_type* _ptr_nonlocal_var;
};





//template <typename T, std::size_t Dim, typename Container>
//class small_strain_isotropic_nonlocal_damage_damage_variable :
//        public history_material_base<T>,
//        public small_strain_material_base<T, Dim, Container>,
//        public nonlocal_material_base<T, Dim>
//{
//public:
//    using value_type = T;
//    using size_type = std::size_t;
//    using tensor4 = tmech::tensor<value_type, Dim, 4>;

//    small_strain_isotropic_nonlocal_damage_damage_variable():
//        history_material_base<T>(2),
//        small_strain_material_base<T, Dim, Container>(),
//        nonlocal_material_base<T, Dim>(1ul),
//        _material(nullptr),
//        _damage_func(nullptr),
//        _yield_function(nullptr)
//    {}

//    small_strain_isotropic_nonlocal_damage_damage_variable(small_strain_material_base<T, Dim, Container> * material, propagation_law_base<T> * damage_func, yield_function_isotropic_damage_strain_based<T, Dim, Container> * yield_function):
//        history_material_base<T>(2),
//        small_strain_material_base<T, Dim, Container>(),
//        nonlocal_material_base<T, Dim>(1ul),
//        _material(material),
//        _damage_func(damage_func),
//        _yield_function(yield_function)
//    {}

//    virtual ~small_strain_isotropic_nonlocal_damage_damage_variable(){}

//    inline void init()override{
//        if(!this->_is_init){
//            _material->init();
//            this->_C = _material->tangent_tensor();
//            this->_is_init = true;
//        }
//    }

//    inline void update()override{
//        //easy handling
//        auto& eps{this->_strain};
//        auto& sig{this->_stress};
//        auto& C{this->_C};
//        auto& history{this->_history};

//        //update base material
//        //compute trail stress
//        _material->strain_tensor() = eps;
//        _material->update();

//        //update equivalent scalar
//        _yield_function->update_equivalent_scalar();

//        //solve for new damage multiplier
//        const auto his_new{_yield_function->solve(history[1])};

//        if(_yield_function->yielding(history[1]) /*&& this->nonlocal_variables[0] > history[0]*/){
//            _yield_function->update_strain_dependent_parts();
//            const auto D{std::max(this->_nonlocal_variables[0], history[0])};//{this->nonlocal_variables[0]};//{std::max(this->nonlocal_variables[0], history[0])};
//            const auto dD{_damage_func->derivative(_yield_function->state())};
//            C = (1-D)*_material->tangent_tensor();
//            sig = (1-D)*_material->stress_tensor();
//            this->_receiver[0] = _material->stress_tensor();
//            this->_source[0] = dD*_yield_function->derivative();
//            history[0] = D;
//            history[1] = his_new;
//            this->_element_marker = true;
//        }else{
//            //No damage
//            const auto D{history[0]};
//            C = (1-D)*_material->tangent_tensor();
//            sig = (1-D)*_material->stress_tensor();
//            this->_receiver[0].fill(0);
//            this->_source[0].fill(0);
//            this->_element_marker = false;
//            history[1] = std::max(_yield_function->critical_value(), history[1]);
//        }
//    }

//    inline void update_stress()override{
//        throw std::runtime_error("small_strain_isotropic_nonlocal_d_damage::update_stress()");
//    }

//    inline void update_tangent()override{
//        throw std::runtime_error("small_strain_isotropic_nonlocal_d_damage::update_tangent()");
//    }

//    constexpr inline auto set_base_material(small_strain_material_base<value_type, Dim, Container> * __material){
//        _material = __material;
//    }

//    constexpr inline auto set_propagation_law(propagation_law_base<value_type> * __propagation_law){
//        _damage_func = __propagation_law;
//    }

//    constexpr inline auto set_yield_function(yield_function_isotropic_damage_strain_based<value_type, Dim, Container> * __yield_function){
//        _yield_function = __yield_function;
//    }


//    virtual inline material_base<T, Dim, Container>* base_material(){
//        return _material->base_material();
//    }

//    inline void update_nonlocal_variables() override {
//        _material->strain_tensor() = this->_strain;
//        _yield_function->update_equivalent_scalar();
//        if(_yield_function->yielding(this->_history[1])){
//            this->_nonlocal_variables[0] = _damage_func->value(_yield_function->state());
//        }else {
//            this->_nonlocal_variables[0] = this->_history[0];
//        }
//    }

//protected:
//    small_strain_material_base<value_type, Dim, Container> * _material;
//    propagation_law_base<value_type> * _damage_func;
//    yield_function_isotropic_damage_strain_based<value_type, Dim, Container> * _yield_function;
//};


#endif // SMALL_STRAIN_ISOTROPIC_NONLOCAL_DAMAGE_BONES_H
