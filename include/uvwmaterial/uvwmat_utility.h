/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef UVWMAT_UTILITY_H
#define UVWMAT_UTILITY_H

template <typename T>
struct is_map_container
{
    static constexpr bool value = false;
};

template <typename Key, typename Value>
struct is_map_container<std::map<Key, Value>>
{
    static constexpr bool value = true;
};

template <typename T>
struct is_vector_any_container
{
    static constexpr bool value = false;
};

template <>
struct is_vector_any_container<std::vector<std::any>>
{
    static constexpr bool value = true;
};

template <typename T>
struct is_vector_value_type_container
{
    static constexpr bool value = false;
};

template <typename T>
struct is_vector_value_type_container<std::vector<T>>
{
    static constexpr bool value = true;
};

template <typename Container>
struct set_parameter_in_container
{
    template<typename ...Parameter>
    static constexpr inline auto set_parameter(Container & container, Parameter && ...parameter){
        container.reserve(sizeof...(parameter));
        push_back_parameter(container, parameter...);
    }

private:
    template<typename T, typename ...Parameter>
    static constexpr inline auto push_back_parameter(Container & container, T && first, Parameter && ...parameter){
        container.push_back(first);
        push_back_parameter(container, parameter...);
    }

    static constexpr inline auto push_back_parameter(Container & /*container*/){}
};


template <typename Key, typename Value>
struct set_parameter_in_container<std::map<Key, Value>>
{
    template<typename ...Parameter>
    static constexpr inline auto set_parameter(std::map<Key, Value> & container, Parameter && ...parameter){
        push_back_parameter(container, parameter...);
    }

private:
    template<typename T, typename ...Parameter>
    static constexpr inline auto push_back_parameter(std::map<Key, Value> & container, T && first, Parameter && ...parameter){
        container.insert(first);
        push_back_parameter(container, parameter...);
    }

    static constexpr inline auto push_back_parameter(std::map<Key, Value> & /*container*/){}
};

#endif // UVWMAT_UTILITY_H
