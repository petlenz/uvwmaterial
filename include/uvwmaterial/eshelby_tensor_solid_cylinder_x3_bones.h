#ifndef ESHELBY_TENSOR_SOLID_CYLINDER_X3_BONES_H
#define ESHELBY_TENSOR_SOLID_CYLINDER_X3_BONES_H

template <typename _T, std::size_t _Dim>
class eshelby_tensor_solid_cylinder_x3 : public eshelby_tensor_solid_cylinder<_T, _Dim>
{
public:
    using value_type = _T;
    using base_class = eshelby_tensor_solid_cylinder<_T, _Dim>;

    eshelby_tensor_solid_cylinder_x3();

    virtual~eshelby_tensor_solid_cylinder_x3();
};

#endif // ESHELBY_TENSOR_SOLID_CYLINDER_X3_BONES_H
