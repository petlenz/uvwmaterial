#ifndef SMALL_STRAIN_SHORT_FIBRE_NONLOCAL_COMPOSITE_BONES_H
#define SMALL_STRAIN_SHORT_FIBRE_NONLOCAL_COMPOSITE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class small_strain_short_fibre_damage_nonlocal_composite :
        public small_strain_short_fibre_composite<_T, _Dim, _Container>,
        public short_fibre_history_base<_T, _Dim, _Container>,
        public short_fibre_nonlocal_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr small_strain_short_fibre_damage_nonlocal_composite():
        small_strain_short_fibre_composite<_T, _Dim, _Container>(),
        short_fibre_history_base<_T, _Dim, _Container>(this),
        short_fibre_nonlocal_base<_T, _Dim, _Container>(this)
    {}

    virtual ~small_strain_short_fibre_damage_nonlocal_composite(){}

    inline virtual void init()override{
        small_strain_short_fibre_composite<_T, _Dim, _Container>::init();
        short_fibre_history_base<_T, _Dim, _Container>::init();
        short_fibre_nonlocal_base<_T, _Dim, _Container>::init();
    }

    inline virtual void reinit()override{
        small_strain_short_fibre_composite<_T, _Dim, _Container>::reinit();
    }

    inline virtual void update()override{
        this->set_local_history();
        this->set_internal_nonlocal_variables();
        small_strain_short_fibre_composite<_T, _Dim, _Container>::update();
        short_fibre_nonlocal_base<_T, _Dim, _Container>::update();
        this->get_local_history();
    }

    inline virtual void update_stress() override{
        this->set_local_history();
        this->set_internal_nonlocal_variables();
        small_strain_short_fibre_composite<_T, _Dim, _Container>::update_stress();
        this->get_local_history();
    }

    inline virtual void update_tangent() override{
        this->set_local_history();
        this->set_internal_nonlocal_variables();
        small_strain_short_fibre_composite<_T, _Dim, _Container>::update_tangent();
        this->get_local_history();
    }
};

#endif // SMALL_STRAIN_SHORT_FIBRE_NONLOCAL_COMPOSITE_BONES_H
