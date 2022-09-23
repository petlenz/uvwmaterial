#ifndef LINEAR_ISOTROPIC_THERMAL_CONDUCTIVITY_BONES_H
#define LINEAR_ISOTROPIC_THERMAL_CONDUCTIVITY_BONES_H

//#include <tmech/tmech.h>
#include "thermal_material_base_bones.h"

template <typename _T, std::size_t _Dim, typename _Container>
class linear_isotropic_thermal_conductivity : public conductivity_material_base<_T, _Dim, _Container>
{
public:
    linear_isotropic_thermal_conductivity() {}

    virtual ~linear_isotropic_thermal_conductivity() {}

    inline void init(){
        if(!this->_is_init){
            const tmech::eye<_T, _Dim, 2> I;
            this->_conductivity = this->_parameter[0]*I;
            this->_is_init = true;
        }
    }

    inline void reinit(){
        this->_is_init = false;
        this->init();
    }

    inline void update(){
        const tmech::eye<_T, _Dim, 2> I;
        this->_flux = -this->_conductivity*this->_gradient;
        this->_dflux.fill(0);
    }

    inline void update_flux(){
        this->_flux = -this->_conductivity*this->_gradient;
    }

    inline void update_conductivity(){
        //do nothing here
    }

    inline conductivity_material_base<_T, _Dim, _Container>* base_material(){
        return this;
    }

private:
};
#endif // LINEAR_ISOTROPIC_THERMAL_CONDUCTIVITY_BONES_H
