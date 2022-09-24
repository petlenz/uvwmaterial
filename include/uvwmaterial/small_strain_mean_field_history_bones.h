/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SMALL_STRAIN_MEAN_FIELD_HISTORY_BONES_H
#define SMALL_STRAIN_MEAN_FIELD_HISTORY_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class small_strain_mean_field_composite_history :
        public small_strain_mean_field_composite<_T, _Dim, _Container>,
        public composite_history_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr small_strain_mean_field_composite_history():
        small_strain_mean_field_composite<_T, _Dim, _Container>(),
        composite_history_base<_T, _Dim, _Container>(this)
    {}

    virtual ~small_strain_mean_field_composite_history(){}

    inline virtual void init() override {
        small_strain_mean_field_composite<_T, _Dim, _Container>::init();
        composite_history_base<_T, _Dim, _Container>::init();
    }

    inline virtual void reinit() override {
        small_strain_mean_field_composite<_T, _Dim, _Container>::reinit();
    }

    inline virtual void update()override{
        this->set_local_history();
        small_strain_mean_field_composite<_T, _Dim, _Container>::update();
        this->get_local_history();
    }

    inline virtual void update_stress() override{
        this->set_local_history();
        small_strain_mean_field_composite<_T, _Dim, _Container>::update_stress();
        this->get_local_history();
    }

    inline virtual void update_tangent() override{
        this->set_local_history();
        small_strain_mean_field_composite<_T, _Dim, _Container>::update_tangent();
        this->get_local_history();
    }
};

#endif // SMALL_STRAIN_MEAN_FIELD_HISTORY_BONES_H
