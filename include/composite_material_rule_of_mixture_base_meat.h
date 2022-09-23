#ifndef COMPOSITE_MATERIAL_RULE_OF_MIXTURE_BASE_MEAT_H
#define COMPOSITE_MATERIAL_RULE_OF_MIXTURE_BASE_MEAT_H

template<typename _T, std::size_t _Dim, typename _Container>
constexpr rule_of_mixture_composite_base<_T, _Dim, _Container>::rule_of_mixture_composite_base():
    _materials(),
    _volume_fractions(),
    _kernal(nullptr)
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr rule_of_mixture_composite_base<_T, _Dim, _Container>::rule_of_mixture_composite_base(rule_of_mixture_kernal_base<_T, _Dim, _Container> * __kernal):
    _materials(),
    _volume_fractions(),
    _kernal(__kernal)
{}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr rule_of_mixture_composite_base<_T, _Dim, _Container>::rule_of_mixture_composite_base(rule_of_mixture_kernal_base<_T, _Dim, _Container> *  __kernal, size_type const __number_of_history):
    _materials(__number_of_history),
    _volume_fractions(__number_of_history),
    _kernal(__kernal)
{}

template<typename _T, std::size_t _Dim, typename _Container>
rule_of_mixture_composite_base<_T, _Dim, _Container>::~rule_of_mixture_composite_base(){}

template<typename _T, std::size_t _Dim, typename _Container>
inline void rule_of_mixture_composite_base<_T, _Dim, _Container>::init(){
    if(_kernal == nullptr){
        throw std::runtime_error("small_strain_rule_of_mixture_base: no kernal");
    }

    if(_materials.empty()){
        throw std::runtime_error("small_strain_rule_of_mixture_base: no materials");
    }

    for(auto material : _materials){
        material->init();
    }

    this->_strain_concentration_tensors.resize(this->_materials.size());
    _kernal->determine_strain_concentration_tensors(_materials,_volume_fractions,this->_strain_concentration_tensors);
    composite_material_base<_T, _Dim, _Container>::init_materials();
}

template<typename _T, std::size_t _Dim, typename _Container>
inline void rule_of_mixture_composite_base<_T, _Dim, _Container>::reinit(){
    for(auto material : _materials){
        material->init();
    }
    _kernal->determine_strain_concentration_tensors(_materials,_volume_fractions,this->_strain_concentration_tensors);
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto rule_of_mixture_composite_base<_T, _Dim, _Container>::push_back(small_strain_material_base<_T, _Dim, _Container> * __material){
    _materials.push_back(__material);
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto rule_of_mixture_composite_base<_T, _Dim, _Container>::push_back(value_type const __volume_fraction){
    _volume_fractions.push_back(__volume_fraction);
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto rule_of_mixture_composite_base<_T, _Dim, _Container>::push_back(value_type const __volume_fraction, small_strain_material_base<_T, _Dim, _Container> * __material){
    _materials.push_back(__material);
    _volume_fractions.push_back(__volume_fraction);
}

template<typename _T, std::size_t _Dim, typename _Container>
inline typename rule_of_mixture_composite_base<_T, _Dim, _Container>::value_type const& rule_of_mixture_composite_base<_T, _Dim, _Container>::volume_fraction(size_type __idx_material)const{
    return _volume_fractions[__idx_material];
}

template<typename _T, std::size_t _Dim, typename _Container>
inline typename rule_of_mixture_composite_base<_T, _Dim, _Container>::value_type& rule_of_mixture_composite_base<_T, _Dim, _Container>::volume_fraction(size_type __idx_material){
    return _volume_fractions[__idx_material];
}

template<typename _T, std::size_t _Dim, typename _Container>
inline typename rule_of_mixture_composite_base<_T, _Dim, _Container>::size_type rule_of_mixture_composite_base<_T, _Dim, _Container>::number_of_materials() const{
    return _materials.size();
}

template<typename _T, std::size_t _Dim, typename _Container>
inline material_base<_T, _Dim, _Container> const* rule_of_mixture_composite_base<_T, _Dim, _Container>::material(size_type __idx_material)const{
    return _materials[__idx_material];
}

template<typename _T, std::size_t _Dim, typename _Container>
inline material_base<_T, _Dim, _Container>* rule_of_mixture_composite_base<_T, _Dim, _Container>::material(size_type __idx_material){
    return _materials[__idx_material];
}

template<typename _T, std::size_t _Dim, typename _Container>
constexpr inline auto rule_of_mixture_composite_base<_T, _Dim, _Container>::set_kernal(rule_of_mixture_kernal_base<_T, _Dim, _Container> * __kernal){
    _kernal = __kernal;
}

#endif // COMPOSITE_MATERIAL_RULE_OF_MIXTURE_BASE_MEAT_H
