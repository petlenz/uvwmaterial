#ifndef INCLUSION_MEAT_H
#define INCLUSION_MEAT_H

template <typename _T, std::size_t _Dim, typename _Container>
inclusion<_T, _Dim, _Container>::inclusion():
    inclusion_base<_T, _Dim, _Container>()
{}

template <typename _T, std::size_t _Dim, typename _Container>
inclusion<_T, _Dim, _Container>::inclusion(value_type const c, material_base<_T, _Dim, _Container> * material, eshelby_tensor_base<_T, _Dim> * S):
    inclusion_base<_T, _Dim, _Container>(c, material, S)
{}

template <typename _T, std::size_t _Dim, typename _Container>
inclusion<_T, _Dim, _Container>::inclusion(inclusion const& data):
    inclusion_base<_T, _Dim, _Container>(data.c, data.material, data.S)
{}

template <typename _T, std::size_t _Dim, typename _Container>
inclusion<_T, _Dim, _Container>::~inclusion(){}

#endif // INCLUSION_MEAT_H
