#ifndef SMALL_STRAIN_MATERIAL_BASE_BONES_H
#define SMALL_STRAIN_MATERIAL_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class small_strain_material_base :
        public solid_material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;

    small_strain_material_base();

    small_strain_material_base(std::vector<value_type> const& __parameter);

    small_strain_material_base(std::initializer_list<value_type> const& __parameter);

    virtual ~small_strain_material_base(){}

    template<typename ..._Parameter>
    small_strain_material_base(_Parameter&& ... __parameter);

    constexpr inline auto& strain_tensor();

    constexpr inline auto const& strain_tensor() const;

    virtual inline void init() = 0;

    virtual inline void reinit() = 0;

    virtual inline void update() = 0;

    virtual inline void update_stress() = 0;

    virtual inline void update_tangent() = 0;

protected:
    tmech::tensor<value_type, _Dim, 2> _strain;
};

#endif // SMALL_STRAIN_MATERIAL_BASE_BONES_H
