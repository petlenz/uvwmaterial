/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef CONDUCTIVITY_MATERIAL_BASE_BONES_H
#define CONDUCTIVITY_MATERIAL_BASE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class conductivity_material_base : public material_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr conductivity_material_base();

    constexpr conductivity_material_base(_Container const& __parameter);

    template<typename _Parameter>
    constexpr conductivity_material_base(std::initializer_list<_Parameter> const& __parameter);

    template<typename ..._Parameter>
    constexpr conductivity_material_base(_Parameter&& ... __parameter);

    virtual ~conductivity_material_base(){}

    constexpr inline auto const& dflux_tensor() const{
        return _dflux;
    }

    constexpr inline auto& dflux_tensor(){
        return _dflux;
    }

    constexpr inline auto const& flux_tensor() const;

    constexpr inline auto& flux_tensor();

    constexpr inline auto const& conductivity_tensor() const;

    constexpr inline auto& conductivity_tensor();

    constexpr inline auto& gradient_tensor();

    constexpr inline auto const& gradient_tensor()const;

    constexpr inline auto const& scalar_value()const;

    constexpr inline auto& scalar_value();

    inline virtual void init() = 0;

    inline virtual void reinit() = 0;

    inline virtual void update() = 0;

    inline virtual void update_flux() = 0;

    inline virtual void update_conductivity() = 0;

    inline virtual conductivity_material_base<_T, _Dim, _Container>* base_material() = 0;

protected:
    value_type                 _scalar;
    tmech::tensor<_T, _Dim, 1> _flux;
    tmech::tensor<_T, _Dim, 1> _dflux;
    tmech::tensor<_T, _Dim, 2> _conductivity;
    tmech::tensor<_T, _Dim, 1> _gradient;
};





#endif // CONDUCTIVITY_MATERIAL_BASE_BONES_H
