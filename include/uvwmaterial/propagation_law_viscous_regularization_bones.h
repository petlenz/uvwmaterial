/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef PROPAGATION_LAW_VISCOUS_REGULARIZATION_BONES_H
#define PROPAGATION_LAW_VISCOUS_REGULARIZATION_BONES_H

template <typename _T>
class propagation_law_viscous_regularization : public propagation_law_base<_T>
{
public:
    using value_type = _T;

    propagation_law_viscous_regularization();

    propagation_law_viscous_regularization(propagation_law_base<_T> * __base_law);

    virtual ~propagation_law_viscous_regularization(){}

    inline value_type value(value_type const __gamma)const override;

    inline value_type derivative(value_type const __gamma)const override;

    constexpr inline auto set_parameter(value_type const __eta);

    constexpr inline auto set_history(value_type const __history);

    constexpr inline auto set_time_increment(value_type const __dtime);

    constexpr inline auto set_base_law(propagation_law_base<_T> & __base_law);

private:
    propagation_law_base<value_type> * _base_law;
    value_type _eta;
    value_type _history;
    value_type _dtime;
};


#endif // PROPAGATION_LAW_VISCOUS_REGULARIZATION_BONES_H
