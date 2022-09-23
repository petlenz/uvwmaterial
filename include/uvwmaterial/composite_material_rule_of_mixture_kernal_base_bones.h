#ifndef COMPOSITE_MATERIAL_RULE_OF_MIXTURE_KERNAL_BASE_BONES_H
#define COMPOSITE_MATERIAL_RULE_OF_MIXTURE_KERNAL_BASE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class rule_of_mixture_kernal_base
{
public:

    rule_of_mixture_kernal_base(){}

    virtual ~rule_of_mixture_kernal_base(){}

    constexpr inline virtual void determine_strain_concentration_tensors(
            std::vector<small_strain_material_base<_T, _Dim, _Container>*> const& _materials,
            std::vector<_T> const& _volume_fractions,
            std::vector<tmech::tensor<_T, _Dim, 4>> & _A) const = 0;
};


//voigt
template <typename _T, std::size_t _Dim, typename _Container>
class rule_of_mixture_voigt : public rule_of_mixture_kernal_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    rule_of_mixture_voigt();

    virtual ~rule_of_mixture_voigt();

    constexpr inline virtual void determine_strain_concentration_tensors(std::vector<small_strain_material_base<_T, _Dim, _Container>*> const& _materials,
                                                                         std::vector<_T> const& _volume_fractions,
                                                                         std::vector<tmech::tensor<_T, _Dim, 4>> & _A) const override;
};

//reuss
template <typename _T, std::size_t _Dim, typename _Container>
class rule_of_mixture_reuss : public rule_of_mixture_kernal_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    rule_of_mixture_reuss();

    virtual ~rule_of_mixture_reuss();

    constexpr inline virtual void determine_strain_concentration_tensors(
            std::vector<small_strain_material_base<_T, _Dim, _Container>*> const& _materials,
            std::vector<_T> const& _volume_fractions,
            std::vector<tmech::tensor<_T, _Dim, 4>> & _A) const override;

};


//vrh
template <typename _T, std::size_t _Dim, typename _Container>
class rule_of_mixture_vrh : public rule_of_mixture_kernal_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    rule_of_mixture_vrh();

    virtual ~rule_of_mixture_vrh();

    constexpr inline virtual void determine_strain_concentration_tensors(std::vector<small_strain_material_base<_T, _Dim, _Container>*> const& _materials,
                                                                         std::vector<_T> const& _volume_fractions,
                                                                         std::vector<tmech::tensor<_T, _Dim, 4>> & _A) const override;
};
#endif // COMPOSITE_MATERIAL_RULE_OF_MIXTURE_KERNAL_BASE_BONES_H
