#ifndef COMPOSITE_MATERIAL_INCLUSION_BASE_BONES_H
#define COMPOSITE_MATERIAL_INCLUSION_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class inclusion_base
{
public:
    using value_type = _T;

    inclusion_base();

    inclusion_base(value_type const __c, material_base<_T, _Dim, _Container> * __material, eshelby_tensor_base<_T, _Dim> * __S);

    virtual ~inclusion_base(){}

    inline void init();

    constexpr inline auto const& volume_fraction()const;

    constexpr inline auto& volume_fraction();

    constexpr inline auto material()const;

    constexpr inline auto material();

    constexpr inline auto material(material_base<value_type, _Dim, _Container> * __material);

    constexpr inline auto eshelby_tensor()const;

    constexpr inline auto eshelby_tensor();

    constexpr inline auto eshelby_tensor(eshelby_tensor_base<value_type, _Dim> * __S);

    constexpr inline bool is_initialized()const;

    constexpr inline auto reinit(){
        _material->reinit();
        _S->reinit();
    }

protected:
    value_type _c;
    material_base<value_type, _Dim, _Container> * _material;
    eshelby_tensor_base<value_type, _Dim> * _S;
    bool _is_init;

private:
    constexpr inline auto check_data();
};


#endif // COMPOSITE_MATERIAL_INCLUSION_BASE_BONES_H
