#ifndef ESHELBY_TENSOR_BASE_MEAT_H
#define ESHELBY_TENSOR_BASE_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_base<_T, _Dim>::eshelby_tensor_base():
    _parameter(),
    _is_init(false)
{}

template <typename _T, std::size_t _Dim>
constexpr inline bool eshelby_tensor_base<_T, _Dim>::is_initialized()const{
    return _is_init;
}

template <typename _T, std::size_t _Dim>
template<typename ...Parameter>
constexpr inline auto eshelby_tensor_base<_T, _Dim>::set_parameter(Parameter && ...__parameter){
    detail::set_parameter_in_container<std::vector<value_type>>::set_parameter(_parameter, __parameter...);
}
#endif // ESHELBY_TENSOR_BASE_MEAT_H
