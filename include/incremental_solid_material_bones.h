#ifndef INCREMENTAL_SOLID_MATERIAL_BONES_H
#define INCREMENTAL_SOLID_MATERIAL_BONES_H

template <typename _T, std::size_t _Dim>
class incremental_solid_material
{
public:
    using value_type = _T;

    incremental_solid_material();

    constexpr inline auto const& dstrain_tensor()const;

    constexpr inline auto& dstrain_tensor();

protected:
    tmech::tensor<_T, _Dim, 2> _dstrain;
};


#endif // INCREMENTAL_SOLID_MATERIAL_BONES_H
