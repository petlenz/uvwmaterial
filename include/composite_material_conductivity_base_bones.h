#ifndef COMPOSITE_MATERIAL_CONDUCTIVITY_BASE_BONES_H
#define COMPOSITE_MATERIAL_CONDUCTIVITY_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class composite_material_conductivity_base :
        public composite_material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    composite_material_conductivity_base();

    virtual ~composite_material_conductivity_base();

    inline virtual value_type& volume_fraction(size_type const __idx) = 0;

    inline virtual value_type const& volume_fraction(size_type const __idx) const = 0;

    inline virtual size_type number_of_materials() const = 0;

    inline virtual material_base<_T, _Dim, _Container> const* material(size_type const __idx) const = 0;

    inline virtual material_base<_T, _Dim, _Container>* material(size_type const __idx) = 0;

    constexpr inline auto const& gradient_concentration_tensors()const;

    constexpr inline auto& gradient_concentration_tensors();

    inline virtual void update_gradient() = 0;    

protected:
    std::vector<tmech::tensor<value_type, _Dim, 2>> _gradient_concentration_tensors;
};

#endif // COMPOSITE_MATERIAL_CONDUCTIVITY_BASE_BONES_H
