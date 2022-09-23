#ifndef EIGENSTRAIN_MATERIAL_BASE_MEAT_H
#define EIGENSTRAIN_MATERIAL_BASE_MEAT_H

template <typename _T, std::size_t _Dim>
eigenstrain_material_base<_T, _Dim>::eigenstrain_material_base():
    _eigenstrain_tensor()
{}

template <typename _T, std::size_t _Dim>
eigenstrain_material_base<_T, _Dim>::~eigenstrain_material_base(){}

template <typename _T, std::size_t _Dim>
constexpr inline auto const& eigenstrain_material_base<_T, _Dim>::eigenstrain_tensor()const{
    return _eigenstrain_tensor;
}

template <typename _T, std::size_t _Dim>
constexpr inline auto& eigenstrain_material_base<_T, _Dim>::eigenstrain_tensor(){
    return _eigenstrain_tensor;
}

#endif // EIGENSTRAIN_MATERIAL_BASE_MEAT_H
