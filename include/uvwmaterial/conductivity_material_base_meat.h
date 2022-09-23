#ifndef CONDUCTIVITY_MATERIAL_BASE_MEAT_H
#define CONDUCTIVITY_MATERIAL_BASE_MEAT_H


template<typename _T, std::size_t _Dim, typename _Container>
constexpr conductivity_material_base<_T, _Dim, _Container>::conductivity_material_base():
    material_base<_T, _Dim, _Container>(),
    _scalar(),
    _flux(),
    _conductivity(),
    _gradient()
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr conductivity_material_base<_T, _Dim, _Container>::conductivity_material_base(_Container const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _scalar(),
    _flux(),
    _conductivity(),
    _gradient()
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename _Parameter>
constexpr conductivity_material_base<_T, _Dim, _Container>::conductivity_material_base(std::initializer_list<_Parameter> const& __parameter):
    material_base<_T, _Dim, _Container>(__parameter),
    _scalar(),
    _flux(),
    _conductivity(),
    _gradient()
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename ..._Parameter>
constexpr conductivity_material_base<_T, _Dim, _Container>::conductivity_material_base(_Parameter&& ...__parameter):
    material_base<_T, _Dim, _Container>(__parameter...),
    _scalar(),
    _flux(),
    _conductivity(),
    _gradient()
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& conductivity_material_base<_T, _Dim, _Container>::flux_tensor()const{
    return _flux;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& conductivity_material_base<_T, _Dim, _Container>::flux_tensor(){
    return _flux;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& conductivity_material_base<_T, _Dim, _Container>::conductivity_tensor()const{
    return _conductivity;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& conductivity_material_base<_T, _Dim, _Container>::conductivity_tensor(){
    return _conductivity;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& conductivity_material_base<_T, _Dim, _Container>::gradient_tensor()const{
    return _gradient;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& conductivity_material_base<_T, _Dim, _Container>::gradient_tensor(){
    return _gradient;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& conductivity_material_base<_T, _Dim, _Container>::scalar_value()const{
    return _scalar;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& conductivity_material_base<_T, _Dim, _Container>::scalar_value(){
    return _scalar;
}
#endif // CONDUCTIVITY_MATERIAL_BASE_MEAT_H
