/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef TIME_DEPENDENT_MATERIAL_BASE_BONES_H
#define TIME_DEPENDENT_MATERIAL_BASE_BONES_H

template <typename T>
class time_dependent_material_base
{
public:
    using value_type = T;

    time_dependent_material_base();

    virtual ~time_dependent_material_base();

    //inline virtual void reinit() = 0;

    constexpr inline auto& time();

    constexpr inline auto const& time()const;

    constexpr inline auto& delta_time();

    constexpr inline auto const& delta_time()const;

protected:
    value_type _time;
    value_type _delta_time;
};

#endif // TIME_DEPENDENT_MATERIAL_BASE_BONES_H
