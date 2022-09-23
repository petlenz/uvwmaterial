#ifndef THERMO_MECHANICAL_MATERIAL_BASE_BONES_H
#define THERMO_MECHANICAL_MATERIAL_BASE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class thermo_mechanical_material_base :
        public material_base<_T, _Dim, _Container>,
        public temperature_dependent_material_base<_T>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr thermo_mechanical_material_base();

    constexpr thermo_mechanical_material_base(_Container const& __parameter);

    template<typename _Parameter>
    constexpr thermo_mechanical_material_base(std::initializer_list<_Parameter> const& __parameter);

    template<typename ..._Parameter>
    constexpr thermo_mechanical_material_base(_Parameter&& ... __parameter);

    virtual ~thermo_mechanical_material_base(){}

    constexpr inline auto const& heat_flux_tensor() const;

    constexpr inline auto& heat_flux_tensor();

    constexpr inline auto const& dflux_tensor() const{
        return _dq;
    }

    constexpr inline auto& dflux_tensor(){
        return _dq;
    }

    constexpr inline auto const& conductivity_tensor() const;

    constexpr inline auto& conductivity_tensor();

    constexpr inline auto const& thermal_gradient() const;

    constexpr inline auto& thermal_gradient();

    constexpr inline auto const& stress_tensor() const;

    constexpr inline auto& stress_tensor();

    constexpr inline auto const& tangent_tensor() const;

    constexpr inline auto& tangent_tensor();

    constexpr inline auto& temperature();

    constexpr inline auto const& temperature()const;

    constexpr inline auto& delta_temperature();

    constexpr inline auto const& delta_temperature()const;

    constexpr inline auto& heat_capacity(){
        return _cd;
    }

    constexpr inline auto const& heat_capacity()const{
        return _cd;
    }

    constexpr inline auto& density(){
        return _rho;
    }

    constexpr inline auto const& density()const{
        return _rho;
    }

    constexpr inline auto& latent_heat(){
        return _latent_heat;
    }

    constexpr inline auto const& latent_heat()const{
        return _latent_heat;
    }

    constexpr inline auto& dlatent_heat(){
        return _dlatent_heat;
    }

    constexpr inline auto const& dlatent_heat()const{
        return _dlatent_heat;
    }

    inline virtual void init() = 0;

    inline virtual void reinit() = 0;

    inline virtual void update() = 0;

    inline virtual void update_stress() = 0;

    inline virtual void update_tangent() = 0;

    inline virtual void update_heat_flux() = 0;

    inline virtual void update_conductivity() = 0;

    inline virtual thermo_mechanical_material_base<_T, _Dim, _Container>* base_material() = 0;

protected:
    tmech::tensor<_T, _Dim, 2> _stress; //elastic stress
    tmech::tensor<_T, _Dim, 4> _C; //material tangent
    //tmech::tensor<_T, _Dim, 2> _eig_strain; //eigenstrain tensor thermal
    tmech::tensor<_T, _Dim, 1> _q; //heat flux
    tmech::tensor<_T, _Dim, 2> _k; //thermal conductivity
    tmech::tensor<_T, _Dim, 1> _dT; //thermal gradient
    tmech::tensor<_T, _Dim, 1> _dq; //dq/dT
    value_type _cd;
    value_type _rho;
    value_type _temperature;
    value_type _delta_temperature;
    value_type _latent_heat;
    value_type _dlatent_heat;
};

#endif // THERMO_MECHANICAL_MATERIAL_BASE_BONES_H
