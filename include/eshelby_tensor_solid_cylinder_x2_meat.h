#ifndef ESHELBY_TENSOR_SOLID_CYLINDER_X2_MEAT_H
#define ESHELBY_TENSOR_SOLID_CYLINDER_X2_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder_x2<_T, _Dim>::eshelby_tensor_solid_cylinder_x2()
{
    base_class::direction() = tmech::tensor<_T, 3, 1>{0,1,0};
}

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_cylinder_x2<_T, _Dim>::~eshelby_tensor_solid_cylinder_x2(){}

#endif // ESHELBY_TENSOR_SOLID_CYLINDER_X2_MEAT_H
