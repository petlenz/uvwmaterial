#ifndef ESHELBY_TENSOR_NONLINEAR_BASE_BONES_H
#define ESHELBY_TENSOR_NONLINEAR_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class eshelby_tensor_nonlinear
{
public:
    eshelby_tensor_nonlinear():
        _base_material(nullptr),
        _isotropization(nullptr),
        _parameter()
    {}

    virtual ~eshelby_tensor_nonlinear(){}

    constexpr inline auto set_base_material(material_base<_T, _Dim, _Container>* __base_material){
        _base_material = __base_material;
    }

    constexpr inline auto set_isotropization(tensor_isotropization_base<_T, _Dim, _Container>* __isotropization){
        _isotropization = __isotropization;
    }

    constexpr inline auto init(){
        _isotropization->init();
    }

//    template<typename ..._Parameter>
//    constexpr inline void set_parameter(_Parameter... __parameter){
//        detail::set_parameter_in_container<_Container>::set_parameter(_parameter, __parameter...);
//    }

protected:
    material_base<_T, _Dim, _Container>* _base_material;
    tensor_isotropization_base<_T, _Dim, _Container>* _isotropization;
    _Container _parameter;
};

#endif // ESHELBY_TENSOR_NONLINEAR_BASE_BONES_H
