#ifndef PROPAGATION_LAW_BASE_BONES_H
#define PROPAGATION_LAW_BASE_BONES_H


template <typename T>
class propagation_law_base
{
public:
    using value_type = T;

    propagation_law_base():parameter(){}

    template<typename ...Parameter>
    propagation_law_base(Parameter && ...__parameter):parameter() {
        detail::set_parameter_in_container<std::vector<value_type>>::set_parameter(parameter, __parameter...);
    }

    propagation_law_base(propagation_law_base const& data):parameter(data.parameter) {}

    virtual ~propagation_law_base(){}

    virtual inline value_type value(value_type const gamma) const = 0;

    virtual inline value_type derivative(value_type const gamma) const = 0;

    template<typename ...Parameter>
    constexpr inline void set_parameter(Parameter &&... __parameter){
        detail::set_parameter_in_container<std::vector<value_type>>::set_parameter(parameter, __parameter...);
    }

protected:
    std::vector<value_type> parameter;
};


#endif // PROPAGATION_LAW_BASE_BONES_H
