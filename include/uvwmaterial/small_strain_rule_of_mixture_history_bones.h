/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SMALL_STRAIN_RULE_OF_MIXTURE_HISTORY_BONES_H
#define SMALL_STRAIN_RULE_OF_MIXTURE_HISTORY_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class small_strain_rule_of_mixture_history_base
        : public rule_of_mixture_history_base<_T, _Dim, _Container>, public history_material_base<_T>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr small_strain_rule_of_mixture_history_base(){}

    virtual ~small_strain_rule_of_mixture_history_base(){}

    inline virtual void update()override{
        this->set_local_history();
        small_strain_rule_of_mixture_base<_T, _Dim, _Container>::update();
        this->get_local_history();
    }

    inline virtual void update_stress() override{
        this->set_local_history();
        small_strain_rule_of_mixture_base<_T, _Dim, _Container>::update_stress();
        this->get_local_history();
    }

    inline virtual void update_tangent() override{
        this->set_local_history();
        small_strain_rule_of_mixture_base<_T, _Dim, _Container>::update_tangent();
        this->get_local_history();
    }

    constexpr inline virtual void init()override{
        small_strain_rule_of_mixture_history_base<_T, _Dim, _Container>::init();
    }
};


#endif // SMALL_STRAIN_RULE_OF_MIXTURE_HISTORY_BONES_H
