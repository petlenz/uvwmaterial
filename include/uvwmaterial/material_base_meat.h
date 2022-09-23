#ifndef MATERIAL_BASE_MEAT_H
#define MATERIAL_BASE_MEAT_H


template<typename _T, std::size_t _Dim, typename _Container>
constexpr material_base<_T, _Dim, _Container>::material_base():
    _parameter(),
    _is_init(false)
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr material_base<_T, _Dim, _Container>::material_base(_Container const& __parameter):
    _parameter(__parameter),
    _is_init(false)
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename _Parameter>
constexpr material_base<_T, _Dim, _Container>::material_base(std::initializer_list<_Parameter> const& __parameter):
    _parameter(__parameter),
    _is_init(false)
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename ..._Parameter>
constexpr material_base<_T, _Dim, _Container>::material_base(_Parameter&& ...__parameter):
    _parameter(),
    _is_init(false)
{
    detail::set_parameter_in_container<_Container>::set_parameter(_parameter, __parameter...);
}


template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline void material_base<_T, _Dim, _Container>::set_parameter(_Container const& __parameter){
    _parameter = __parameter;
}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename ..._Parameter>
constexpr inline void material_base<_T, _Dim, _Container>::set_parameter(_Parameter... __parameter){
    detail::set_parameter_in_container<_Container>::set_parameter(_parameter, __parameter...);
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& material_base<_T, _Dim, _Container>::parameter()const{
    return _parameter;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& material_base<_T, _Dim, _Container>::parameter(){
    return _parameter;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline bool material_base<_T, _Dim, _Container>::is_initialized()const{
    return _is_init;
}

#endif // MATERIAL_BASE_MEAT_H
