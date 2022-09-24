/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef HISTORY_MATERIAL_BASE_BONES_H
#define HISTORY_MATERIAL_BASE_BONES_H

template<typename _T>
class history_material_base
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr history_material_base();

    constexpr history_material_base(size_type const __number_of_history);

    virtual ~history_material_base();

    template<typename IterBegin, typename IterEnd>
    constexpr inline auto set_history(IterBegin begin, IterEnd end);

    constexpr inline auto& history();

    constexpr inline auto const& history()const;

    constexpr inline auto begin(){
        return _history.begin();
    }

    constexpr inline auto end(){
        return _history.end();
    }

    constexpr inline auto size()const;

    constexpr inline auto resize(size_type const __number_of_history);

protected:
    std::vector<value_type> _history;
};

#endif // HISTORY_MATERIAL_BASE_BONES_H
