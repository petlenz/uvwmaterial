#ifndef COMPOSITE_MATERIAL_EIGENSTRAIN_BASE_MEAT_H
#define COMPOSITE_MATERIAL_EIGENSTRAIN_BASE_MEAT_H

template <typename _T, std::size_t _Dim, typename _Container>
composite_eigenstrain_material_base<_T, _Dim, _Container>::composite_eigenstrain_material_base(material_base<_T, _Dim, _Container> * __material):
    _material(__material),
    _stress_concentration_tensors(),
    _eigenstrain_concentration_tensors()
{}

template <typename _T, std::size_t _Dim, typename _Container>
composite_eigenstrain_material_base<_T, _Dim, _Container>::composite_eigenstrain_material_base(material_base<_T, _Dim, _Container> * __material, size_type const __number_of_inclusions, size_type const __number_of_eigenstrains):
    _material(__material),
    _stress_concentration_tensors(__number_of_inclusions+1),
    _eigenstrain_concentration_tensors(__number_of_eigenstrains)
{
    for(auto& eigenstrain : _eigenstrain_concentration_tensors){
        eigenstrain.resize(__number_of_inclusions+1);
    }
}

#endif // COMPOSITE_MATERIAL_EIGENSTRAIN_BASE_MEAT_H
