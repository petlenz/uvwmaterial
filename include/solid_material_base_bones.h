#ifndef SOLID_MATERIAL_BASE_BONES_H
#define SOLID_MATERIAL_BASE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class solid_material_base : public material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr solid_material_base();

    constexpr solid_material_base(_Container const& __parameter);

    template<typename _Parameter>
    constexpr solid_material_base(std::initializer_list<_Parameter> const& __parameter);

    template<typename ..._Parameter>
    constexpr solid_material_base(_Parameter&& ... __parameter);

    virtual ~solid_material_base(){}

    constexpr inline auto const& stress_tensor() const;

    constexpr inline auto& stress_tensor();

    constexpr inline auto const& tangent_tensor() const;

    constexpr inline auto& tangent_tensor();

    inline virtual void init() = 0;

    inline virtual void reinit() = 0;

    inline virtual void update() = 0;

    inline virtual void update_stress() = 0;

    inline virtual void update_tangent() = 0;

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() = 0;

protected:
    tmech::tensor<_T, _Dim, 2> _stress;
    tmech::tensor<_T, _Dim, 4> _C;
};

#endif // SOLID_MATERIAL_BASE_BONES_H
