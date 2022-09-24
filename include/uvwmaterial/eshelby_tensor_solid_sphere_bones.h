/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_SPHERE_BONES_H
#define ESHELBY_TENSOR_SOLID_SPHERE_BONES_H

template <typename T, std::size_t Dim>
class eshelby_tensor_solid_sphere : public eshelby_tensor_solid_base<T, Dim>
{
public:
    using value_type = T;

    eshelby_tensor_solid_sphere();

    virtual ~eshelby_tensor_solid_sphere(){}

    virtual inline void init() override;

    virtual inline void reinit(){}
};


#endif // ESHELBY_TENSOR_SOLID_SPHERE_BONES_H
