#ifndef HISTORY_MATERIAL_BASE_MEAT_H
#define HISTORY_MATERIAL_BASE_MEAT_H

template<typename T>
constexpr history_material_base<T>::history_material_base():
    _history()
{}

template<typename T>
constexpr history_material_base<T>::history_material_base(size_type const __number_of_history):
    _history(__number_of_history, 0)
{}

template<typename T>
history_material_base<T>::~history_material_base(){}

template<typename T>
template<typename IterBegin, typename IterEnd>
constexpr inline auto history_material_base<T>::set_history(IterBegin begin, IterEnd end){
    std::copy(begin, end, _history.begin());
}

template<typename T>
constexpr inline auto const& history_material_base<T>::history()const{
    return _history;
}

template<typename T>
constexpr inline auto& history_material_base<T>::history(){
    return _history;
}

template<typename T>
constexpr inline auto history_material_base<T>::size()const{
    return _history.size();
}

template<typename T>
constexpr inline auto history_material_base<T>::resize(size_type const __number_of_history){
    _history.resize(__number_of_history);
}

#endif // HISTORY_MATERIAL_BASE_MEAT_H
