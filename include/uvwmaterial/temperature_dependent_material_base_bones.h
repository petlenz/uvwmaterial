#ifndef TEMPERATURE_DEPENDENT_MATERIAL_BASE_BONES_H
#define TEMPERATURE_DEPENDENT_MATERIAL_BASE_BONES_H

template <typename T>
class temperature_dependent_material_base
{
public:
    using value_type = T;

    temperature_dependent_material_base();

    virtual ~temperature_dependent_material_base();

    inline virtual void reinit() = 0;

    constexpr inline auto& temperature();

    constexpr inline auto const& temperature()const;

    constexpr inline auto& delta_temperature();

    constexpr inline auto const& delta_temperature()const;

protected:
    value_type _temperature;
    value_type _delta_temperature;
};

#endif // TEMPERATURE_DEPENDENT_MATERIAL_BASE_BONES_H
