#ifndef STATE_FUNCTION_BASE_MEAT_H
#define STATE_FUNCTION_BASE_MEAT_H

template <typename T, std::size_t Dim, typename Container>
state_function_base<T, Dim, Container>::state_function_base():
    _material(nullptr)
{}

template <typename T, std::size_t Dim, typename Container>
state_function_base<T, Dim, Container>::state_function_base(solid_material_base<value_type, Dim, Container> * __material):
    _material(__material)
{}

template <typename T, std::size_t Dim, typename Container>
constexpr inline auto state_function_base<T, Dim, Container>::set_base_material(solid_material_base<value_type, Dim, Container> * __material){
    _material = __material;
}


////------------------------------------------------------------------------------------------------------------
//template <std::size_t Dim, typename Container>
//state_function_base<float, Dim, Container>::state_function_base():
//    material(nullptr)
//{}

//template <std::size_t Dim, typename Container>
//state_function_base<float, Dim, Container>::state_function_base(material_base<value_type, Dim, Container> * _material):
//    material(_material)
//{}

//template <std::size_t Dim, typename Container>
//constexpr inline auto state_function_base<float, Dim, Container>::set_base_material(material_base<value_type, Dim, Container> * __material){
//    material = __material;
//}

#endif // STATE_FUNCTION_BASE_MEAT_H
