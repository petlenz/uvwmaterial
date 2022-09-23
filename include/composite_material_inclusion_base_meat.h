#ifndef COMPOSITE_MATERIAL_INCLUSION_BASE_MEAT_H
#define COMPOSITE_MATERIAL_INCLUSION_BASE_MEAT_H

template <typename _T, std::size_t _Dim, typename _Container>
inclusion_base<_T, _Dim, _Container>::inclusion_base():
    _c(0),
    _material(nullptr),
    _S(nullptr),
    _is_init(false)
{}

template <typename _T, std::size_t _Dim, typename _Container>
inclusion_base<_T, _Dim, _Container>::inclusion_base(value_type const __c, material_base<_T, _Dim, _Container> * __material, eshelby_tensor_base<_T, _Dim> * __S):
    _c(__c),
    _material(__material),
    _S(__S),
    _is_init(false)
{}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& inclusion_base<_T, _Dim, _Container>::volume_fraction()const{
    return _c;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& inclusion_base<_T, _Dim, _Container>::volume_fraction(){
    _is_init = false;
    return _c;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto inclusion_base<_T, _Dim, _Container>::material()const{
    return _material;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto inclusion_base<_T, _Dim, _Container>::material(material_base<value_type, _Dim, _Container> * __material){
    _is_init = false;
    _material = __material;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto inclusion_base<_T, _Dim, _Container>::eshelby_tensor()const{
    return _S;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto inclusion_base<_T, _Dim, _Container>::eshelby_tensor(){
    return _S;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto inclusion_base<_T, _Dim, _Container>::eshelby_tensor(eshelby_tensor_base<value_type, _Dim> * __S){
    _is_init = false;
    _S = __S;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline bool inclusion_base<_T, _Dim, _Container>::is_initialized()const{
    return _is_init;
}

template <typename _T, std::size_t _Dim, typename _Container>
inline void inclusion_base<_T, _Dim, _Container>::init(){
    if(!_is_init){
        check_data();
        _material->init();
        _S->init();
        _is_init = true;
    }
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto inclusion_base<_T, _Dim, _Container>::material(){
    return _material;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto inclusion_base<_T, _Dim, _Container>::check_data(){
    if(!_material){
        throw std::runtime_error("inclusion_base::check_data(): material is not set");
    }
    if(!_S){
        throw std::runtime_error("inclusion_base::check_data(): eshelby tensor is not set");
    }
}

#endif // COMPOSITE_MATERIAL_INCLUSION_BASE_MEAT_H
