#ifndef ESHELBY_TENSOR_SOLID_CYLINDER_BONES_H
#define ESHELBY_TENSOR_SOLID_CYLINDER_BONES_H

template <typename _T, std::size_t _Dim>
class eshelby_tensor_solid_cylinder : public eshelby_tensor_solid_base<_T, _Dim>
{
public:
    using value_type = _T;
    using base_class = eshelby_tensor_solid_base<_T, _Dim>;

    eshelby_tensor_solid_cylinder();

    virtual ~eshelby_tensor_solid_cylinder();

    virtual inline void init() override;

    inline virtual void reinit(){}

    constexpr inline auto& direction();

    constexpr inline auto const& direction()const;

private:
    template<typename _Tensor, typename _N>
    static constexpr inline auto setup_entries(_Tensor& _S, _T const __nue, _N const& __n);

    tmech::tensor<value_type, 3, 1> _n;
};

#endif // ESHELBY_TENSOR_SOLID_CYLINDER_BONES_H
