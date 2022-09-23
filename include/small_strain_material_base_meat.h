#ifndef SMALL_STRAIN_MATERIAL_BASE_MEAT_H
#define SMALL_STRAIN_MATERIAL_BASE_MEAT_H

template<typename _T, std::size_t _Dim, typename _Container>
small_strain_material_base<_T, _Dim, _Container>::small_strain_material_base():
    _strain()
{}

template<typename _T, std::size_t _Dim, typename _Container>
small_strain_material_base<_T, _Dim, _Container>::small_strain_material_base(std::vector<value_type> const& __parameter):
    solid_material_base<_T, _Dim, _Container>(__parameter),
    _strain()
{}

template<typename _T, std::size_t _Dim, typename _Container>
small_strain_material_base<_T, _Dim, _Container>::small_strain_material_base(std::initializer_list<value_type> const& __parameter):
    solid_material_base<_T, _Dim, _Container>(__parameter),
    _strain()
{}

template<typename _T, std::size_t _Dim, typename _Container>
template<typename ..._Parameter>
small_strain_material_base<_T, _Dim, _Container>::small_strain_material_base(_Parameter&& ... __parameter):
    solid_material_base<_T, _Dim, _Container>(__parameter...),
    _strain()
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& small_strain_material_base<_T, _Dim, _Container>::strain_tensor(){
    return _strain;
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& small_strain_material_base<_T, _Dim, _Container>::strain_tensor()const{
    return _strain;
}

#endif // SMALL_STRAIN_MATERIAL_BASE_MEAT_H
