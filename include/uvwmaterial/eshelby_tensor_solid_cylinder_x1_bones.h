#ifndef ESHELBY_TENSOR_SOLID_CYLINDER_X1_BONES_H
#define ESHELBY_TENSOR_SOLID_CYLINDER_X1_BONES_H

template <typename _T, std::size_t _Dim>
class eshelby_tensor_solid_cylinder_x1 :
        public eshelby_tensor_solid_cylinder<_T, _Dim>
{
public:
    using value_type = _T;
    using base_class = eshelby_tensor_solid_cylinder<_T, _Dim>;

    eshelby_tensor_solid_cylinder_x1();

    virtual~eshelby_tensor_solid_cylinder_x1();
};

#endif // ESHELBY_TENSOR_SOLID_CYLINDER_X1_BONES_H
