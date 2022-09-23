#ifndef THERMAL_MATERIAL_BASE_BONES_H
#define THERMAL_MATERIAL_BASE_BONES_H

#include "material_base_bones.h"

template<typename _T, std::size_t _Dim, typename _Container>
class thermal_material_base : public material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr thermal_material_base();

    constexpr thermal_material_base(_Container const& __parameter);

    template<typename _Parameter>
    constexpr thermal_material_base(std::initializer_list<_Parameter> const& __parameter);

    template<typename ..._Parameter>
    constexpr thermal_material_base(_Parameter&& ... __parameter);

    virtual ~thermal_material_base(){}

    constexpr inline auto const& heat_flux_tensor() const;

    constexpr inline auto& heat_flux_tensor();

    constexpr inline auto const& conductivity_tensor() const;

    constexpr inline auto& conductivity_tensor();

    constexpr inline auto const& thermal_gradient() const;

    constexpr inline auto& thermal_gradient();

    inline virtual void init() = 0;

    inline virtual void update() = 0;

    inline virtual void update_heat_flux() = 0;

    inline virtual void update_conductivity() = 0;

    inline virtual thermal_material_base<_T, _Dim, _Container>* base_material() = 0;

protected:
    tmech::tensor<_T, _Dim, 1> _q; //heat flux
    tmech::tensor<_T, _Dim, 2> _k; //thermal conductivity
    tmech::tensor<_T, _Dim, 1> _dT; //thermal gradient
};


#endif // THERMAL_MATERIAL_BASE_BONES_H
