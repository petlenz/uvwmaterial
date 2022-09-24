/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef MATERIAL_BASE_BONES_H
#define MATERIAL_BASE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class material_base
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr material_base();

    constexpr material_base(_Container const& __parameter);

    template<typename _Parameter>
    constexpr material_base(std::initializer_list<_Parameter> const& __parameter);

    template<typename ..._Parameter>
    constexpr material_base(_Parameter&& ... __parameter);

    virtual ~material_base(){}

    constexpr inline void set_parameter(_Container const& __parameter);

    template<typename ..._Parameter>
    constexpr inline void set_parameter(_Parameter... __parameter);

    constexpr inline auto const& parameter() const;

    constexpr inline auto& parameter();

    constexpr inline bool is_initialized()const;

    inline virtual void init() = 0;

    inline virtual void reinit() = 0;

protected:
    _Container _parameter;
    bool _is_init;
};

#endif // MATERIAL_BASE_BONES_H
