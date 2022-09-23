#ifndef ESHELBY_TENSOR_SOLID_SPHERE_MEAT_H
#define ESHELBY_TENSOR_SOLID_SPHERE_MEAT_H

#include "eshelby_tensor_functions.h"

template <typename T, std::size_t Dim>
eshelby_tensor_solid_sphere<T, Dim>::eshelby_tensor_solid_sphere()
{}

template <typename T, std::size_t Dim>
inline void eshelby_tensor_solid_sphere<T, Dim>::init() {
    //Micromechanics-based homogenization of the effective
    //physical properties of composites with an anisotropic
    //matrix and interfacial imperfections
    //const auto I{tmech::eye<double, 3, 2>()};
    //S1 = (5*nue-1)/(15*(1-nue))*tmech::otimes(I, I) + (4-5*nue)/(15*(1-nue)) *(boxtimes(I,I) + transr(transi(otimes(I,I))));

    if(!this->_is_init){
        if constexpr (Dim == 2){
            tmech::tensor<value_type, 3, 4> Stemp;
            eshelby_tensor_sphere_func(Stemp, this->_parameter[0]);
            this->_S = convert_3D_to_2D(Stemp);
        }else{
            eshelby_tensor_sphere_func(this->_S, this->_parameter[0]);
        }
        this->_is_init = true;
    }
}

#endif // ESHELBY_TENSOR_SOLID_SPHERE_MEAT_H
