#ifndef MEAN_FIELD_COMPOSITE_KERNAL_BASE_BONES_H
#define MEAN_FIELD_COMPOSITE_KERNAL_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_base;

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_kernal_base
{
public:
    mean_field_composite_kernal_base() {}

    virtual ~mean_field_composite_kernal_base(){}

    constexpr inline auto set_material(mean_field_composite_base<_T, _Dim, _Container>* __material){
        _material = __material;
    }
protected:
    mean_field_composite_base<_T, _Dim, _Container>* _material;
};

#endif // MEAN_FIELD_COMPOSITE_KERNAL_BASE_BONES_H
