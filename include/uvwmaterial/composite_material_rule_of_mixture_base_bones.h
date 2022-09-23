#ifndef COMPOSITE_MATERIAL_RULE_OF_MIXTURE_BASE_BONES_H
#define COMPOSITE_MATERIAL_RULE_OF_MIXTURE_BASE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class rule_of_mixture_composite_base :
        public composite_material_solid_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr rule_of_mixture_composite_base();

    constexpr rule_of_mixture_composite_base(rule_of_mixture_kernal_base<_T, _Dim, _Container> * __kernal);

    constexpr rule_of_mixture_composite_base(rule_of_mixture_kernal_base<_T, _Dim, _Container> *  __kernal, size_type const __number_of_history);

    virtual ~rule_of_mixture_composite_base();

    inline void init();

    inline void reinit();

    constexpr inline auto push_back(small_strain_material_base<_T, _Dim, _Container> * __material);

    constexpr inline auto push_back(value_type const __volume_fraction);

    constexpr inline auto push_back(value_type const __volume_fraction, small_strain_material_base<_T, _Dim, _Container> * __material);

    inline virtual value_type const& volume_fraction(size_type __idx_material)const;

    inline virtual value_type& volume_fraction(size_type __idx_material);

    inline virtual size_type number_of_materials() const override;

    inline virtual material_base<_T, _Dim, _Container> const* material(size_type __idx_material)const override;

    inline virtual material_base<_T, _Dim, _Container>* material(size_type __idx_material) override;

    constexpr inline auto set_kernal(rule_of_mixture_kernal_base<_T, _Dim, _Container> * __kernal);

protected:
    std::vector<small_strain_material_base<_T, _Dim, _Container>*> _materials;
    std::vector<value_type> _volume_fractions;
private:
    rule_of_mixture_kernal_base<_T, _Dim, _Container> * _kernal;
};

#endif // COMPOSITE_MATERIAL_RULE_OF_MIXTURE_BASE_BONES_H
