#ifndef INCLUSION_BONES_H
#define INCLUSION_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class inclusion :
        public inclusion_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    inclusion();

    inclusion(value_type const c, material_base<_T, _Dim, _Container> * material, eshelby_tensor_base<_T, _Dim> * S);

    inclusion(inclusion const& data);

    virtual ~inclusion();
};

#endif // INCLUSION_BONES_H
