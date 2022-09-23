#ifndef ESHELBY_TENSOR_CONDUCTIVITY_BASE_MEAT_H
#define ESHELBY_TENSOR_CONDUCTIVITY_BASE_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_conductivity_base<_T, _Dim>::eshelby_tensor_conductivity_base():
    _S()
{}

template <typename _T, std::size_t _Dim>
eshelby_tensor_conductivity_base<_T, _Dim>::~eshelby_tensor_conductivity_base(){}

template <typename _T, std::size_t _Dim>
inline constexpr auto const& eshelby_tensor_conductivity_base<_T, _Dim>::tensor()const{
    return _S;
}

#endif // ESHELBY_TENSOR_CONDUCTIVITY_BASE_MEAT_H
