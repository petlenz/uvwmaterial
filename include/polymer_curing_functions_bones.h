// Copyright 2021 Peter Lenz
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------
#ifndef POLYMER_CURING_FUNCTIONS_BONES_H
#define POLYMER_CURING_FUNCTIONS_BONES_H

//Kinetic models describing the chemical kinetics of thermoset cure
//Karkanas, Partridge, Attwood 1996
//Modelling the cure of a commercial epoxy resin for applications in resin transfer moulding

template<typename T>
class kinetic_model_thermoset_base
{
public:
    using value_type = T;

    kinetic_model_thermoset_base() = delete ;

    template<typename Derived>
    kinetic_model_thermoset_base(Derived & derived_class):
        delta_z(0),
        z_n(0),
        time_inc(0),
        theta(0),
        parameter(),
        derived(&derived_class)
    {}

    template<typename Derived, typename ...Parameter>
    kinetic_model_thermoset_base(Derived & derived_class, Parameter ... __parameter):
        delta_z(0),
        z_n(0),
        time_inc(0),
        theta(0),
        parameter(),
        derived(&derived_class)
    {
        parameter.reserve(sizeof...(__parameter));
        push_back_parameter(__parameter...);
    }

    virtual ~kinetic_model_thermoset_base(){}


    constexpr inline std::pair<T, T> operator()(value_type const& __delta_z){
        delta_z = __delta_z;
        return (*this)();
        //std::make_tuple(derived->value(), derived->derivative());
    }

    template<typename Parameter>
    constexpr inline void set_parameter(Parameter const& __parameter){
        push_back_parameter(__parameter);
    }

    template<typename ...Parameter>
    constexpr inline void set_parameter(Parameter ... __parameter){
        parameter.reserve(sizeof...(__parameter));
        push_back_parameter(__parameter...);
    }

    constexpr inline auto set_history(value_type const __history){
        z_n = __history;
    }

    constexpr inline auto set_time_increment(value_type const __time_inc){
        time_inc = __time_inc;
    }

    constexpr inline auto set_temperature(value_type const __theta){
        theta = __theta;
    }

    constexpr inline auto get_degree_of_cure()const{
        return z_n + delta_z;
    }

    constexpr inline auto update(){
        delta_z = 0;
        if(z_n < 1.0){
            //std::cout<<"z_n = "<<z_n<<std::endl;
            delta_z = boost::math::tools::newton_raphson_iterate(std::ref(*this), 5e-5, 0., 1., std::numeric_limits<value_type>::digits);
            //std::cout<<"z_n = "<<z_n<<" dz = "<<delta_z<<std::endl;
        }
    };

    virtual inline value_type value() = 0;

    virtual inline value_type derivative() = 0;

    virtual inline std::pair<T, T> operator()() const = 0;

private:
    template<typename ...Parameter>
    constexpr inline auto push_back_parameter(value_type const First, Parameter ...__parameter){
        parameter.push_back(First);
        push_back_parameter(__parameter...);
    }

    constexpr inline auto push_back_parameter(value_type const Last){
        parameter.push_back(Last);
    }

protected:
    value_type delta_z;
    value_type z_n;
    value_type time_inc;
    value_type theta;
    std::vector<value_type> parameter;
    kinetic_model_thermoset_base<value_type> * derived;
    static constexpr value_type R{8.31446261815324};
};



template <typename T>
class nth_order_reaction : public kinetic_model_thermoset_base<T>
{
public:
    using value_type = T;
    using pair = std::pair<value_type, value_type>;

    nth_order_reaction(value_type const A, value_type const E, value_type const n):
        kinetic_model_thermoset_base<T>(*this, A, E, n)
    {}

    inline auto operator()(value_type const curing, value_type const theta)const{
        const auto A{this->parameter[0]};
        const auto E{this->parameter[1]};
        const auto n{this->parameter[2]};
        const value_type k{A*std::exp(-E/(this->R*theta))};
        return k*std::pow((1-curing),n);
    }

    inline pair operator()(){
        const auto A{this->parameter[0]};
        const auto E{this->parameter[1]};
        const auto n{this->parameter[2]};
        const value_type k{A*std::exp(-E/(this->R*this->theta))};
        //r = dz - dt*K*std::pow((1-(z+dz)),n)
        const auto r{this->delta_z - this->time_inc*k*std::pow((1-(this->z_n + this->delta_z)), n)};
        const auto dr{1. + k*n*this->time_inc*pow(1-(this->z_n + this->delta_z),n-1)};
        return std::make_tuple(r, dr);
    }
};


template <typename T>
class autocatalytic_reaction : public kinetic_model_thermoset_base<T>
{
public:
    using value_type = T;
    using pair = std::pair<value_type, value_type>;

    autocatalytic_reaction():
        kinetic_model_thermoset_base<T>(*this)
    {}

    autocatalytic_reaction(value_type const A, value_type const E, value_type const n, value_type const m):
        kinetic_model_thermoset_base<T>(*this, A, E, n, m)
    {}

    inline auto operator()(value_type const curing, value_type const theta)const{
        const auto A{this->parameter[0]};
        const auto E{this->parameter[1]};
        const auto n{this->parameter[2]};
        const auto m{this->parameter[3]};

        const value_type k{A*std::exp(-E/(this->R*theta))};
        const value_type nth{k*std::pow((1 - curing),n)};
        return nth*std::pow(curing, m);
    }


    virtual inline pair operator()()const{
        const auto A{this->parameter[0]};
        const auto E{this->parameter[1]};
        const auto n{this->parameter[2]};
        const auto m{this->parameter[3]};

        //r = dz - dt*k*(1 - (z + dz))^n*(z+dz)^m
        const value_type k{A*std::exp(-E/(this->R*this->theta))};
        const value_type r{this->delta_z - this->time_inc*k*std::pow((1 - (this->delta_z + this->z_n)),n)*std::pow(this->delta_z + this->z_n, m)};
        const value_type dr{1.0 + k*this->time_inc*std::pow(this->delta_z + this->z_n, m-1)*std::pow(1-(this->z_n+this->delta_z), n-1)
                    *(m*(this->z_n+this->delta_z-1)+n*(this->z_n+this->delta_z))};
        //std::cout<<r<<" "<<dr<<std::endl;
        return std::pair(r, dr);
    }

    inline value_type value() override {
        const auto A{this->parameter[0]};
        const auto E{this->parameter[1]};
        const auto n{this->parameter[2]};
        const auto m{this->parameter[3]};

        const value_type k{A*std::exp(-E/(this->R*this->theta))};
        const auto z{this->delta_z + this->z_n};
        return k*std::pow(z, m)*std::pow(1.-z, n);
        //return this->delta_z - this->time_inc*k*std::pow((1 - (this->delta_z + this->z_n)),n)*std::pow(this->delta_z + this->z_n, m);

        //auto f = function();
        //f.update(parameter_pack());
        //return f(parameter_pack());
    }

    inline value_type derivative() override {
        const auto A{this->parameter[0]};
        const auto E{this->parameter[1]};
        const auto n{this->parameter[2]};
        const auto m{this->parameter[3]};

        //r = dz - dt*k*(1 - (z + dz))^n*(z+dz)^m
        const value_type k{A*std::exp(-E/(this->R*this->theta))};
        return 1 + k*this->time_inc*std::pow(this->delta_z + this->z_n, m-1)*std::pow(1-(this->z_n+this->delta_z), n-1)
                *(m*(this->z_n+this->delta_z-1)+n*(this->z_n+this->delta_z));
        //return symdiff::derivative<1>(function(), dz)(parameter_pack());
    }

//    constexpr inline auto function(){
//        const auto k{A*symdiff::exp(-E/(c_R*c_theta))};
//        const auto nth{k*symdiff::pow((one - (zn + dz)),n)};
//        return nth*symdiff::pow(zn + dz, m);
//    }

private:
//    constexpr inline auto parameter_pack()const{
//        return std::make_tuple(std::ref(this->delta_z),
//                               std::ref(this->parameter[0]),
//                               std::ref(this->parameter[1]),
//                               std::ref(this->parameter[2]),
//                               std::ref(this->parameter[3]),
//                               std::ref(this->theta),
//                               std::ref(this->R),
//                               std::ref(this->history),
//                               std::ref(this->time_inc));
//    }

//    symdiff::variable<value_type, 0> dz;
//    symdiff::constant<value_type, 1> A;
//    symdiff::constant<value_type, 2> E;
//    symdiff::constant<value_type, 3> m;
//    symdiff::constant<value_type, 4> n;
//    symdiff::constant<value_type, 5> c_theta;
//    symdiff::constant<value_type, 6> c_R;
//    symdiff::constant<value_type, 7> zn;
//    symdiff::constant<value_type, 8> dt;
//    symdiff::real<value_type, 1, 0, 1> one;
};




template <typename T>
class nth_order_autocatalytic_reaction : public kinetic_model_thermoset_base<T>
{
public:
    using value_type = T;
    using pair = std::pair<value_type, value_type>;

    nth_order_autocatalytic_reaction(value_type A1, value_type A2, value_type E1, value_type E2, value_type m, value_type n):
        kinetic_model_thermoset_base<T>(*this, A1, A2, E1, E2, m, n)
    {}

    nth_order_autocatalytic_reaction():
        kinetic_model_thermoset_base<T>(*this)
    {}


    inline pair operator()()const{
        const auto A{this->parameter[0]};
        const auto E{this->parameter[1]};
        const auto n{this->parameter[2]};
        const auto m{this->parameter[3]};

        //r = dz - dt*k*(1 - (z + dz))^n*(z+dz)^m
        const value_type k{A*std::exp(-E/(this->R*this->theta))};
        const value_type r{this->delta_z - this->time_inc*k*std::pow((1 - (this->delta_z + this->z_n)),n)*std::pow(this->delta_z + this->z_n, m)};
        const value_type dr{1 + k*this->time_inc*std::pow(this->delta_z + this->z_n, m-1)*std::pow(1-(this->z_n+this->delta_z), n-1)
                    *(m*(this->z_n+this->delta_z-1)+n*(this->z_n+this->delta_z))};
        return std::pair(r, dr);
    }


//    inline value_type value() override {
//        auto f = function();
//        const auto para{parameter_pack()};
//        f.reset();
//        f.update(para);
//        return f(para);
//    }

//    inline value_type derivative() override {
//        return symdiff::derivative<1>(function(), dz)(parameter_pack());
//    }

//    constexpr inline auto function(){
//        auto k1{A1*symdiff::exp(-E1/(c_R*c_theta))};
//        auto k2{A2*symdiff::exp(-E2/(c_R*c_theta))};
//        return dz - dt*((k1+k2*symdiff::pow(zn + dz, m))*symdiff::pow((one - (zn + dz)), n));
//    }

private:
//    constexpr inline auto parameter_pack()const{
//        return std::make_tuple(std::ref(this->delta_z),
//                               std::ref(this->parameter[0]),
//                               std::ref(this->parameter[1]),
//                               std::ref(this->parameter[2]),
//                               std::ref(this->parameter[3]),
//                               std::ref(this->parameter[4]),
//                               std::ref(this->parameter[5]),
//                               std::ref(this->theta),
//                               std::ref(this->R),
//                               std::ref(this->history),
//                               std::ref(this->time_inc));
//    }

//    symdiff::variable<value_type, 0> dz;
//    symdiff::constant<value_type, 1> A1;
//    symdiff::constant<value_type, 2> A2;
//    symdiff::constant<value_type, 3> E1;
//    symdiff::constant<value_type, 4> E2;
//    symdiff::constant<value_type, 5> m;
//    symdiff::constant<value_type, 6> n;
//    symdiff::constant<value_type, 7> c_theta;
//    symdiff::constant<value_type, 8> c_R;
//    symdiff::constant<value_type, 9> zn;
//    symdiff::constant<value_type, 10> dt;
//    symdiff::real<value_type, 1, 0, 1> one;
};

#endif // POLYMER_CURING_FUNCTIONS_BONES_H
