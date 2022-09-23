#ifndef PROPAGATION_LAW_VISCOUS_REGULARIZATION_MEAT_H
#define PROPAGATION_LAW_VISCOUS_REGULARIZATION_MEAT_H

template <typename _T>
propagation_law_viscous_regularization<_T>::propagation_law_viscous_regularization():
    propagation_law_base<_T>(),
    _base_law(),
    _eta(0)
{}

template <typename _T>
propagation_law_viscous_regularization<_T>::propagation_law_viscous_regularization(propagation_law_base<_T> * __base_law):
    propagation_law_base<_T>(),
    _base_law(&__base_law),
    _eta(0)
{}

template <typename _T>
inline typename propagation_law_viscous_regularization<_T>::value_type propagation_law_viscous_regularization<_T>::value(value_type const __gamma)const{
    return ((_eta/(_eta+_dtime))*_history + ((_dtime/(_eta+_dtime))*_base_law->value(__gamma)));
}

template <typename _T>
inline typename propagation_law_viscous_regularization<_T>::value_type propagation_law_viscous_regularization<_T>::derivative(value_type const __gamma)const{
    return (_dtime/(_eta+_dtime))*_base_law->derivative(__gamma);
}

template <typename _T>
constexpr inline auto propagation_law_viscous_regularization<_T>::set_parameter(value_type const __eta){
    _eta = __eta;
}

template <typename _T>
constexpr inline auto propagation_law_viscous_regularization<_T>::set_history(value_type const __history){
    _history = __history;
}

template <typename _T>
constexpr inline auto propagation_law_viscous_regularization<_T>::set_time_increment(value_type const __dtime){
    _dtime = __dtime;
}

template <typename _T>
constexpr inline auto propagation_law_viscous_regularization<_T>::set_base_law(propagation_law_base<_T> & __base_law){
    _base_law = &__base_law;
}


#endif // PROPAGATION_LAW_VISCOUS_REGULARIZATION_MEAT_H
