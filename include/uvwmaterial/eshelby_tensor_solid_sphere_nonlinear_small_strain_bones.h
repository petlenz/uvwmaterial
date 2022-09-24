/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_SPHERE_NONLINEAR_SMALL_STRAIN_BONES_H
#define ESHELBY_TENSOR_SOLID_SPHERE_NONLINEAR_SMALL_STRAIN_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class eshelby_tensor_solid_sphere_nonlinear_small_strain :
        public eshelby_tensor_solid_sphere<_T, _Dim>,
        public eshelby_tensor_nonlinear<_T, _Dim, _Container>
{
public:
    using eshelby_nonlinear_base = eshelby_tensor_nonlinear<_T, _Dim, _Container>;
    using eshelby_base           = eshelby_tensor_solid_sphere<_T, _Dim>;

    eshelby_tensor_solid_sphere_nonlinear_small_strain(){}

    virtual ~eshelby_tensor_solid_sphere_nonlinear_small_strain(){}

    virtual void init() override {
        eshelby_base::init();
        eshelby_nonlinear_base::init();
    }

    virtual inline void reinit(){}

    virtual void update() override {
        this->_isotropization->update();
        eshelby_base::_parameter[0] = this->_isotropization->data()[0];
        eshelby_base::_is_init = false;
        eshelby_base::init();
    }
};


#endif // ESHELBY_TENSOR_SOLID_SPHERE_NONLINEAR_SMALL_STRAIN_BONES_H
