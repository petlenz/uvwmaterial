/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef YIELD_FUNCTION_BASE_BONES_H
#define YIELD_FUNCTION_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class yield_function_base
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    yield_function_base():
        _history_n(),
        _scalar_eq(),
        _critical_val(),
        _parameter(),
        _equivalent_func(),
        _base_material(nullptr)
    {}

    yield_function_base(size_type const __number_of_state_functions):
        _history_n(__number_of_state_functions),
        _scalar_eq(__number_of_state_functions),
        _critical_val(__number_of_state_functions),
        _parameter(),
        _equivalent_func(__number_of_state_functions),
        _base_material(nullptr)
    {}

    ~yield_function_base(){}

    constexpr inline auto reserve(size_type __size){
        _history_n.reserve(__size);
        _scalar_eq.reserve(__size);
        _critical_val.reserve(__size);
        _equivalent_func.reserve(__size);
    }

    constexpr inline auto resize(size_type __size){
        _history_n.resize(__size);
        _scalar_eq.resize(__size);
        _critical_val.resize(__size);
        _equivalent_func.resize(__size);
    }

    constexpr inline solid_material_base<_T, _Dim, _Container> * base_material()const{
        return _base_material;
    }

    constexpr inline void base_material(solid_material_base<_T, _Dim, _Container> * __base_material){
        _base_material = __base_material;
    }

    constexpr inline state_function_base<_T, _Dim, _Container>* state_function(size_type __idx)const{
        return _equivalent_func[__idx];
    }

    constexpr inline void state_function(size_type __idx, state_function_base<_T, _Dim, _Container>* __state_function){
        _equivalent_func[__idx] = __state_function;
    }

    constexpr inline void push_back(state_function_base<_T, _Dim, _Container>* __state_function){
        _equivalent_func.push_back(__state_function);
    }

    constexpr inline void critical_value(size_type __idx, value_type const __critical_value){
        _critical_val[__idx] = __critical_value;
    }

    constexpr inline auto critical_value(size_type __idx){
        return _critical_val[__idx];
    }

    template<typename ...Parameter>
    constexpr inline auto set_parameter(Parameter && ... __parameter){
        detail::set_parameter_in_container<std::vector<value_type>>::set_parameter(_parameter, __parameter...);
        std::cout<<_parameter[0]<<std::endl;
    }

    virtual inline value_type solve(value_type const& __history) = 0;

    constexpr inline void update_equivalent_scalar(){
        assert(!_equivalent_func.empty());
        for(size_type i{0}; i<_equivalent_func.size(); ++i){
            _equivalent_func[i]->update();
            _scalar_eq[i] = _equivalent_func[i]->value();
        }
    }

    virtual constexpr inline bool yielding(value_type const& __history) const = 0;

protected:
    std::vector<value_type> _history_n;
    std::vector<value_type> _scalar_eq;
    std::vector<value_type> _critical_val;
    std::vector<value_type> _parameter;
    std::vector<state_function_base<_T, _Dim, _Container>*> _equivalent_func;
    solid_material_base<_T, _Dim, _Container> * _base_material;
};



//template <typename _T, std::size_t _Dim, typename _Container>
//class yield_function_isotropic_base :
//        public yield_function_base<_T, _Dim, _Container>
//{
//public:
//    yield_function_isotropic_base():
//        yield_function_base<_T, _Dim, _Container>(1)at
//    {}

//    constexpr inline state_function_base<_T, _Dim, _Container>* state_function()const{
//        return this->_equivalent_func[0];
//    }

//    constexpr inline void state_function(state_function_base<_T, _Dim, _Container>* __state_function){
//        return this->_equivalent_func[0] = __state_function;
//    }

//protected:

//private:

//};


template <typename _T, std::size_t _Dim, typename _Container>
class yield_function_isotropic_damage_strain_based :
        public yield_function_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    yield_function_isotropic_damage_strain_based():
        yield_function_base<_T, _Dim, _Container>(1)
    {}

    virtual ~yield_function_isotropic_damage_strain_based(){}

    virtual inline value_type solve(value_type const& __history) override {
        check_data();
        return std::max(this->_scalar_eq[0], std::max(this->_critical_val[0], __history));
    }

    constexpr inline bool yielding(value_type const& __history)const override {
        check_data();
        return (this->_scalar_eq[0] > std::max(__history, this->_critical_val[0]) ? true : false);
    }

    constexpr inline auto const& derivative()const{
        return _deps_eq;
    }

    constexpr inline auto update_strain_dependent_parts(){
        check_data();
        _deps_eq = this->_equivalent_func[0]->derivative();
    }

    constexpr inline auto critical_value()const{
        return this->_critical_val[0];
    }

    constexpr inline auto& critical_value(){
        return this->_critical_val[0];
    }

    constexpr inline auto state()const{
        return this->_scalar_eq[0];
    }

protected:
    tmech::tensor<value_type, _Dim, 2> _deps_eq;
private:
    constexpr inline auto check_data()const{
        assert(this->_equivalent_func[0] != nullptr);
        assert(this->_base_material != nullptr);
        assert(this->_critical_val[0] != 0);
    }
};



template <typename _T, std::size_t _Dim, typename _Container>
class yield_function_stress_base : public yield_function_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;

    yield_function_stress_base() = delete;

    yield_function_stress_base(yield_function_stress_base<_T, _Dim, _Container> * __base_function):
        yield_function_base<_T, _Dim, _Container>(),
        base_function(__base_function)
    {}

    yield_function_stress_base(yield_function_stress_base<_T, _Dim, _Container> * __base_function,
                               solid_material_base<_T, _Dim, _Container> * __base_material,
                               propagation_law_base<_T> * __hardening_stress):
        yield_function_base<_T, _Dim, _Container>(__base_material),
        hardening_stress(__hardening_stress),
        base_function(__base_function)
    {}

    template<typename ...Parameter>
    yield_function_stress_base(yield_function_stress_base<_T, _Dim, _Container> * __base_function,
                               material_base<_T, _Dim, _Container> * __base_material,
                               propagation_law_base<_T> * __hardening_stress, Parameter && ...__parameter):
        yield_function_base<_T, _Dim, _Container>(__base_material, __parameter...),
        hardening_stress(__hardening_stress),
        base_function(__base_function)
    {}

    virtual ~yield_function_stress_base(){}

    constexpr inline auto set_hardening_stress(propagation_law_base<_T> * __hardening_stress){
        hardening_stress = __hardening_stress;
    }

    virtual constexpr inline value_type value(value_type const ddlambda) const = 0;

    virtual constexpr inline value_type derivative(value_type const ddlambda) const = 0;

    inline value_type solve(value_type const& __history) override {
        this->_history_n[0] = __history;
        delta_lambda = 0;
        try {
            const auto start{(this->_history_n[0] == 0 ? 1e-8 : 0)};
            delta_lambda = boost::math::tools::newton_raphson_iterate(std::ref(*this), start, 0., this->_history_n[0]+1., std::numeric_limits<value_type>::digits);
        }  catch (const std::runtime_error& error) {
            throw error;
        }
        return delta_lambda;
    }

    constexpr inline auto operator()(value_type const ddelta_lambda){
        //            std::cout<<"x:  "<<ddelta_lambda<<std::endl;
        //            std::cout<<"F:  "<<base_function->value(ddelta_lambda)<<std::endl;
        //            std::cout<<"dF: "<<base_function->derivative(ddelta_lambda)<<std::endl;
        //            std::cout<<"his "<<this->history_n<<std::endl;
        return std::make_tuple(base_function->value(ddelta_lambda), base_function->derivative(ddelta_lambda));
    }

    virtual constexpr inline void update_stress_dependent_parts() = 0;

    constexpr inline auto const& derivative_scalar_wrt_strain()const{
        return ddlambda_deps;
    }

    constexpr inline auto const& derivative_normal_wrt_strain()const{
        return dN_deps;
    }

    constexpr inline auto const& normal()const{
        return N;
    }

protected:
    value_type delta_lambda;
    tmech::tensor<value_type, _Dim, 2> N;
    tmech::tensor<value_type, _Dim, 4> dN_deps;
    tmech::tensor<value_type, _Dim, 2> ddlambda_deps;
    propagation_law_base<value_type> * hardening_stress;

private:
    yield_function_stress_base<_T, _Dim, _Container> const* base_function;
};







template <typename _T, std::size_t _Dim, typename _Container>
class j2_yield_function : public yield_function_stress_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;

    j2_yield_function():
        yield_function_stress_base<_T, _Dim, _Container>(this)
    {}

    j2_yield_function(material_base<_T, _Dim, _Container> * __base_material, propagation_law_base<_T> * __hardening_stress):
        yield_function_stress_base<_T, _Dim, _Container>(this, __base_material, __hardening_stress)
    {}

    template<typename ...Parameter>
    j2_yield_function(solid_material_base<_T, _Dim, _Container> * __base_material,
                      propagation_law_base<_T> * __hardening_stress,
                      Parameter && ...__parameter):
        yield_function_stress_base<_T, _Dim, _Container>(this, __base_material, __hardening_stress, __parameter...)
    {}


    constexpr inline bool yielding(value_type const& __history)const override {
        return (this->_scalar_eq[0]  - (this->_critical_val[0] + this->hardening_stress->value(__history)) > 0);
    }

    constexpr inline value_type value(value_type const ddelta_lambda)const override{
        return this->_scalar_eq[0] - 3*this->_parameter[0]*ddelta_lambda - (this->_critical_val[0] + this->hardening_stress->value(this->_history_n[0] + ddelta_lambda));
    }

    constexpr inline value_type derivative(value_type const ddelta_lambda)const override{
        return - 3*this->_parameter[0] - this->hardening_stress->derivative(this->_history_n[0] + ddelta_lambda);
    }

    constexpr inline void update_stress_dependent_parts() override {
        constexpr value_type fac{1.5};
        const auto G{this->_parameter[0]};
        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
        const auto IIvol{tmech::otimes(I, I)/3.};
        const auto IIdev{IIsym - IIvol};
        auto sig_dev{tmech::dev(this->_base_material->stress_tensor())};
        this->N = fac*sig_dev/this->_scalar_eq[0];
        this->dN_deps = G*(3*IIdev - 2*tmech::otimes(this->N, this->N))/this->_scalar_eq[0];
        this->ddlambda_deps = -2*G*this->N/derivative(this->delta_lambda);
    }
};



#endif // YIELD_FUNCTION_BASE_BONES_H
