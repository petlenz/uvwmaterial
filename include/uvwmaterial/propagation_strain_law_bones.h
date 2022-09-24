/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef PROPAGATION_STRAIN_LAW_BONES_H
#define PROPAGATION_STRAIN_LAW_BONES_H


template <typename T>
class damage_exponential1 : public propagation_law_base<T>
{
public:
    using value_type = T;

    damage_exponential1(){}

    damage_exponential1(value_type const __crit_value, value_type const __b):propagation_law_base<T>(__crit_value),b(__b){}

    damage_exponential1(damage_exponential1 const& data):propagation_law_base<T>(data.critical_value),b(data.b){}

    virtual ~damage_exponential1(){}

    inline value_type value(value_type const __gamma)const override{
        return (1.0-std::exp(b*(this->parameter[0] - __gamma)));
    }

     inline value_type derivative(value_type const __gamma)const override{
        return b*std::exp(b*(this->parameter[0] - __gamma));
    }

     constexpr inline auto set_parameter(value_type const __b){
         b = __b;
     }
private:
    value_type b;
};


template <typename T>
class damage_exponential2 : public propagation_law_base<T>
{
public:
    using value_type = T;

    damage_exponential2(){}

    damage_exponential2(value_type const __eps_0, value_type const __eps_f, value_type const __a):
        propagation_law_base<T>(__eps_0),eps_0(__eps_0),eps_f(__eps_f),a(__a){}

    virtual ~damage_exponential2(){}

    inline value_type value(value_type const __gamma)const override{//0.99999
        return a*(static_cast<value_type>(1.0) - (eps_0/__gamma)*std::exp(-(__gamma - eps_0)/(eps_f - eps_0))); //0.99999*
        //return 0.999*(1.0 - std::exp(-(__gamma - eps_0)/(eps_f - eps_0)));
    }

     inline value_type derivative(value_type const __gamma)const override{
         //std::cout<<eps_f<<std::endl;
        //return 0.999*std::exp(-(__gamma - eps_0)/(eps_f - eps_0))/(eps_f - eps_0);
        return a*eps_0*std::exp(-(__gamma - eps_0)/(eps_f - eps_0))*(eps_0 - eps_f - __gamma)/(__gamma*__gamma*(eps_0 - eps_f));
    }

     constexpr inline auto set_parameter(std::vector<value_type> const& __parameter){
         eps_f = __parameter[0];
         eps_0 = __parameter[1];
     }

     constexpr inline auto set_parameter(value_type const __epsf, value_type const __eps0){
         eps_f = __epsf;
         eps_0 = __eps0;
     }
private:
    value_type eps_0;
    value_type eps_f;
    value_type a;
};


template <typename T>
class damage_exponential3 : public propagation_law_base<T>
{
public:
    using value_type = T;

    damage_exponential3(){}

    damage_exponential3(value_type const __kappa0, value_type const __beta, value_type const __alpha):
        propagation_law_base<T>(__kappa0),_kappa0(__kappa0),_beta(__beta),_alpha(__alpha){}

    virtual ~damage_exponential3(){}

    inline value_type value(value_type const __gamma)const override{
        return (1.0 - _kappa0/__gamma*(1.0 - _alpha + _alpha*std::exp(_beta*(_kappa0 - __gamma))))*0.99;
    }

     inline value_type derivative(value_type const __gamma)const override{
         return ((_alpha*_kappa0*((_beta*__gamma+1.0)*std::exp(_beta*(_kappa0-__gamma))-1.0)+_kappa0)/(__gamma*__gamma))*0.99;
    }

private:
    value_type _kappa0;
    value_type _beta;
    value_type _alpha;
};

template <typename T>
class linear : public propagation_law_base<T>
{
public:
    using value_type = T;

    linear():
        K(0)
    {}

    linear(value_type const __K):
        K(__K)
    {}

    constexpr inline value_type value(value_type const __x)const override{
        return K*__x;
    }

    constexpr inline value_type derivative(value_type const /*__x*/)const override{
        return K;
    }

    constexpr inline void set_parameter(value_type const __K){
        K = __K;
    }

private:
    value_type K;
};


template <typename T>
class power_law : public propagation_law_base<T>
{
public:
    using value_type = T;

    power_law():
        K(0),
        m(0)
    {}

    power_law(value_type const __K, value_type const __m):
        K(__K),
        m(__m)
    {}

    constexpr inline value_type value(value_type const __x)const override{
        return K*std::pow(__x, m);
    }

    constexpr inline value_type derivative(value_type const __x)const override{
        return K*m*std::pow(__x, m-1.0);
    }

    constexpr inline void set_parameter(value_type const __K, value_type const __m){
        K = __K;
        m = __m;
    }

private:
    value_type K;
    value_type m;
};


template <typename T>
class exponential_law : public propagation_law_base<T>
{
public:
    using value_type = T;

    exponential_law():
        K(0),
        m(0)
    {}

    exponential_law(value_type const __K, value_type const __m):
        K(__K),
        m(__m)
    {}

    constexpr inline value_type value(value_type const __x)const override{
        return K*(1.0 - std::exp(-m*__x));
    }

    constexpr inline value_type derivative(value_type const __x)const override{
        return K*m*std::exp(-m*__x);
    }

    constexpr inline void set_parameter(value_type const __K, value_type const __m){
        K = __K;
        m = __m;
    }

private:
    value_type K;
    value_type m;
};


#endif // PROPAGATION_STRAIN_LAW_BONES_H
