/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef STATE_FUNCTION_STRESS_BASED_BONES_H
#define STATE_FUNCTION_STRESS_BASED_BONES_H


template<typename T, std::size_t Dim, typename Container>
class von_mises_state_function : public state_function_base<T, Dim, Container>
{
public:
    using value_type = T;

    von_mises_state_function();

    von_mises_state_function(solid_material_base<T, Dim, Container> & material);

    virtual ~von_mises_state_function(){}

    inline void update() override;

    inline value_type value() const override;

    inline tmech::tensor<value_type, Dim, 2> derivative() const override;

};


#endif // STATE_FUNCTION_STRESS_BASED_BONES_H
