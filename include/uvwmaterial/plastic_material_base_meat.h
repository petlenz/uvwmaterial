#ifndef PLASTIC_MATERIAL_BASE_MEAT_H
#define PLASTIC_MATERIAL_BASE_MEAT_H

template <typename _T, std::size_t _Dim>
plastic_material_base<_T, _Dim>::plastic_material_base():
    _inelastic_strain()
{}

template <typename _T, std::size_t _Dim>
plastic_material_base<_T, _Dim>::~plastic_material_base(){}

template <typename _T, std::size_t _Dim>
constexpr inline auto& plastic_material_base<_T, _Dim>::inelastic_strain_tensor(){
    return _inelastic_strain;
}

template <typename _T, std::size_t _Dim>
constexpr inline auto const& plastic_material_base<_T, _Dim>::inelastic_strain_tensor()const{
    return _inelastic_strain;
}

#endif // PLASTIC_MATERIAL_BASE_MEAT_H
