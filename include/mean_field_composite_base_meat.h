#ifndef MEAN_FIELD_COMPOSITE_BASE_MEAT_H
#define MEAN_FIELD_COMPOSITE_BASE_MEAT_H

template <typename _T, std::size_t _Dim, typename _Container>
mean_field_composite_base<_T, _Dim, _Container>::mean_field_composite_base():
    _inclusions(),
    _matrix_material(nullptr),
    _cm(0),
    _kernal(nullptr)
{}

template <typename _T, std::size_t _Dim, typename _Container>
mean_field_composite_base<_T, _Dim, _Container>::~mean_field_composite_base()
{}

template <typename _T, std::size_t _Dim, typename _Container>
inline _T& mean_field_composite_base<_T, _Dim, _Container>::volume_fraction(size_type __idx){
    if(__idx == _inclusions.size()){
        return _cm;
    }else{
        return _inclusions[__idx]->volume_fraction();
    }
}

template <typename _T, std::size_t _Dim, typename _Container>
inline const _T& mean_field_composite_base<_T, _Dim, _Container>::volume_fraction(size_type __idx) const{
    if(__idx == _inclusions.size()){
        return _cm;
    }else{
        return _inclusions[__idx]->volume_fraction();
    }
}

template <typename _T, std::size_t _Dim, typename _Container>
inline typename mean_field_composite_base<_T, _Dim, _Container>::size_type mean_field_composite_base<_T, _Dim, _Container>::number_of_materials() const{
    return _inclusions.size()+1;
}

template <typename _T, std::size_t _Dim, typename _Container>
inline material_base<_T, _Dim, _Container> const* mean_field_composite_base<_T, _Dim, _Container>::material(size_type __idx) const{
    if(__idx == _inclusions.size()){
        return _matrix_material;
    }else{
        return _inclusions[__idx]->material();
    }
}

template <typename _T, std::size_t _Dim, typename _Container>
inline material_base<_T, _Dim, _Container>* mean_field_composite_base<_T, _Dim, _Container>::material(size_type __idx){
    if(__idx == _inclusions.size()){
        return _matrix_material;
    }else{
        return _inclusions[__idx]->material();
    }
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const* mean_field_composite_base<_T, _Dim, _Container>::kernal()const{
    return _kernal;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto mean_field_composite_base<_T, _Dim, _Container>::kernal(mean_field_composite_kernal_base<_T, _Dim, _Container>* __kernal){
    _kernal = __kernal;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const* mean_field_composite_base<_T, _Dim, _Container>::matrix_material()const{
    return _matrix_material;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto mean_field_composite_base<_T, _Dim, _Container>::matrix_material(material_base<_T, _Dim, _Container> * __matrix_material){
    _matrix_material = __matrix_material;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& mean_field_composite_base<_T, _Dim, _Container>::matrix_volume_fraction()const{
    return _cm;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto& mean_field_composite_base<_T, _Dim, _Container>::matrix_volume_fraction(){
    return _cm;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const& mean_field_composite_base<_T, _Dim, _Container>::inclusions()const{
    return _inclusions;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto const* mean_field_composite_base<_T, _Dim, _Container>::inclusion(size_type const __idx)const{
    return _inclusions[__idx];
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto mean_field_composite_base<_T, _Dim, _Container>::inclusion(size_type const __idx, inclusion_base<_T, _Dim, _Container> * __inclusion){
    _inclusions[__idx] = __inclusion;
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto mean_field_composite_base<_T, _Dim, _Container>::push_back(inclusion_base<_T, _Dim, _Container> * __inclusion){
    _inclusions.push_back(__inclusion);
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto mean_field_composite_base<_T, _Dim, _Container>::resize(size_type const __number_of_inclusions){
    _inclusions.resize(__number_of_inclusions);
    this->_strain_concentration_tensors.resize(__number_of_inclusions+1);
}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto mean_field_composite_base<_T, _Dim, _Container>::reserve(size_type const __number_of_inclusions){
    _inclusions.reserve(__number_of_inclusions);
    this->_strain_concentration_tensors.reserve(__number_of_inclusions+1);
}



template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto mean_field_composite_base<_T, _Dim, _Container>::check_data(){
    if(_kernal == nullptr){
        throw std::runtime_error("mean_field_composite_base: no kernal");
    }

    if(_matrix_material == nullptr){
        throw std::runtime_error("mean_field_composite_base: no matrix material");
    }

    if(_cm == 0){
        throw std::runtime_error("mean_field_composite_base: no matrix material volume fraction");
    }

    if(_inclusions.empty()){
        throw std::runtime_error("mean_field_composite_base: no inclusions");
    }
}

#endif // MEAN_FIELD_COMPOSITE_BASE_MEAT_H
