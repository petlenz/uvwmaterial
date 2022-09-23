#ifndef INCREMENTAL_SOLID_MATERIAL_MEAT_H
#define INCREMENTAL_SOLID_MATERIAL_MEAT_H


template <typename _T, std::size_t _Dim>
incremental_solid_material<_T, _Dim>::incremental_solid_material():
    _dstrain()
{}

template <typename _T, std::size_t _Dim>
constexpr inline auto const& incremental_solid_material<_T, _Dim>::dstrain_tensor()const{
    return _dstrain;
}

template <typename _T, std::size_t _Dim>
constexpr inline auto& incremental_solid_material<_T, _Dim>::dstrain_tensor(){
    return _dstrain;
}

#endif // INCREMENTAL_SOLID_MATERIAL_MEAT_H
