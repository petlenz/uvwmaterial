#ifndef EIGENSTRAIN_MATERIAL_BASE_BONES_H
#define EIGENSTRAIN_MATERIAL_BASE_BONES_H

template <typename _T, std::size_t _Dim>
class eigenstrain_material_base
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    eigenstrain_material_base();

    virtual ~eigenstrain_material_base();

    constexpr inline auto const& eigenstrain_tensor()const;

    constexpr inline auto& eigenstrain_tensor();

protected:
    tmech::tensor<value_type, _Dim, 2> _eigenstrain_tensor;
};

#endif // EIGENSTRAIN_MATERIAL_BASE_BONES_H
