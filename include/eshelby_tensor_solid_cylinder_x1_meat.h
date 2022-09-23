#ifndef ESHELBY_TENSOR_SOLID_CYLINDER_X1_MEAT_H
#define ESHELBY_TENSOR_SOLID_CYLINDER_X1_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder_x1<_T, _Dim>::eshelby_tensor_solid_cylinder_x1()
{
    base_class::direction() = tmech::tensor<_T, 3, 1>{1,0,0};
}

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder_x1<_T, _Dim>::~eshelby_tensor_solid_cylinder_x1(){}

#endif // ESHELBY_TENSOR_SOLID_CYLINDER_X1_MEAT_H
