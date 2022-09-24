/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef PLASTIC_MATERIAL_BASE_BONES_H
#define PLASTIC_MATERIAL_BASE_BONES_H

template <typename _T, std::size_t _Dim>
class plastic_material_base : public history_material_base<_T>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    plastic_material_base();

    virtual ~plastic_material_base();

    constexpr inline auto& inelastic_strain_tensor();

    constexpr inline auto const& inelastic_strain_tensor() const;

protected:
    tmech::tensor<value_type, _Dim, 2> _inelastic_strain;
};

#endif // PLASTIC_MATERIAL_BASE_BONES_H
