/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ELEMENT_DELETION_BONES_H
#define ELEMENT_DELETION_BONES_H

template <typename T, std::size_t Dim, typename Container>
class element_deletion
{
public:
    using size_type = std::size_t;

    element_deletion():
        material(nullptr),
        history_material(nullptr),
        history_entry(),
        conditions()
    {}

    constexpr inline auto set_material(material_base<T, Dim, Container>* __material_base){
        material = __material_base;
        history_material = dynamic_cast<history_material_base<T> *>(__material_base);
        assert(history_material);
    }

    constexpr inline auto set_history_entries(std::vector<size_type> const& __history_entry){
        history_entry = __history_entry;
    }

    constexpr inline auto set_conditions(std::vector<T> const& __conditions){
        conditions = __conditions;
    }

    constexpr inline auto set_history_entries(std::vector<size_type> && __history_entry){
        history_entry = std::move(__history_entry);
    }

    constexpr inline auto set_conditions(std::vector<T> && __conditions){
        conditions = std::move(__conditions);
    }

    constexpr inline auto operator()(){
        const auto& his{history_material->history()};
        for(size_type i{0}; i<history_entry.size(); ++i){
            if(his[history_entry[i]] < conditions[i]){
                return false;
            }
        }
        return true;
    }

protected:
    material_base<T, Dim, Container>* material;
    history_material_base<T> * history_material;
    std::vector<size_type> history_entry;
    std::vector<T> conditions;
};

#endif // ELEMENT_DELETION_BONES_H
