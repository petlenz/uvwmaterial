#ifndef MEAN_FIELD_COMPOSITE_BASE_BONES_H
#define MEAN_FIELD_COMPOSITE_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_base //: public composite_material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    mean_field_composite_base();

    virtual ~mean_field_composite_base();

    virtual inline void init() = 0;

    virtual inline void reinit() = 0;

    inline _T& volume_fraction(size_type __idx);

    inline const _T& volume_fraction(size_type __idx) const;

    inline size_type number_of_materials() const;

    inline material_base<_T, _Dim, _Container> const* material(size_type __idx) const;

    inline material_base<_T, _Dim, _Container>* material(size_type __idx);

    constexpr inline auto const* kernal() const;

    constexpr inline auto kernal(mean_field_composite_kernal_base<_T, _Dim, _Container>* __kernal);

    constexpr inline auto const* matrix_material() const;

    constexpr inline auto matrix_material(material_base<_T, _Dim, _Container> * __matrix_material);

    constexpr inline auto const& matrix_volume_fraction() const;

    constexpr inline auto& matrix_volume_fraction();

    constexpr inline auto const& inclusions() const;

    constexpr inline auto const* inclusion(size_type const __idx) const;

    constexpr inline auto inclusion(size_type const __idx, inclusion_base<_T, _Dim, _Container> * __inclusion);

    constexpr inline auto push_back(inclusion_base<_T, _Dim, _Container> * __inclusion);

    constexpr inline auto resize(size_type const __number_of_inclusions);

    constexpr inline auto reserve(size_type const __number_of_inclusions);

protected:
    constexpr inline auto check_data();

    std::vector<inclusion_base<value_type, _Dim, _Container>*> _inclusions;
    material_base<_T, _Dim, _Container> * _matrix_material;
    value_type _cm;
    mean_field_composite_kernal_base<_T, _Dim, _Container> * _kernal;
};

#endif // MEAN_FIELD_COMPOSITE_BASE_BONES_H
