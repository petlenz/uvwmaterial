#ifndef NONLOCAL_MATERIAL_BASE_MEAT_H
#define NONLOCAL_MATERIAL_BASE_MEAT_H

template <typename T, std::size_t Dim>
nonlocal_material_base<T, Dim>::nonlocal_material_base():
    _element_marker(),
    _nonlocal_variables(),
    _source(),
    _receiver()
{}

//template <typename T, std::size_t Dim>
//nonlocal_material_base<T, Dim>::nonlocal_material_base(size_type const __number_of_nonlocal_variables):
//    _element_marker(false),
//    _nonlocal_variables(__number_of_nonlocal_variables),
//    _source(__number_of_nonlocal_variables),
//    _receiver(__number_of_nonlocal_variables)
//{}

//template <typename T, std::size_t Dim>
//nonlocal_material_base<T, Dim>::nonlocal_material_base(std::vector<value_type> const& __nonlocal_variables):
//    _element_marker(false),
//    _nonlocal_variables(__nonlocal_variables),
//    _source(__nonlocal_variables.size()),
//    _receiver(__nonlocal_variables.size())
//{}

template <typename T, std::size_t Dim>
constexpr inline auto const& nonlocal_material_base<T, Dim>::nonlocal_variables()const{
    return _nonlocal_variables;
}

template <typename T, std::size_t Dim>
constexpr inline auto & nonlocal_material_base<T, Dim>::nonlocal_variables(){
    return _nonlocal_variables;
}

#endif // NONLOCAL_MATERIAL_BASE_MEAT_H
