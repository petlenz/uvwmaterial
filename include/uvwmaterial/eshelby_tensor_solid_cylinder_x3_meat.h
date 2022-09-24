/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_CYLINDER_X3_MEAT_H
#define ESHELBY_TENSOR_SOLID_CYLINDER_X3_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder_x3<_T, _Dim>::eshelby_tensor_solid_cylinder_x3()
{
    base_class::direction() = tmech::tensor<_T, 3, 1>{0,0,1};
}

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder_x3<_T, _Dim>::~eshelby_tensor_solid_cylinder_x3(){}

#endif // ESHELBY_TENSOR_SOLID_CYLINDER_X3_MEAT_H
