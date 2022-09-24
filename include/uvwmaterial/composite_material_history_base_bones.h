/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef COMPOSITE_MATERIAL_HISTORY_BASE_BONES_H
#define COMPOSITE_MATERIAL_HISTORY_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class composite_history_base :
        public history_material_base<_T>
{
public:
    using size_type  = std::size_t;
    using value_type = _T;

    composite_history_base(composite_material_base<_T, _Dim, _Container> * __material);

    virtual ~composite_history_base();

    constexpr inline void init();

    inline void set_local_history();

    inline void get_local_history();
private:
    composite_material_base<_T, _Dim, _Container> * _material;
    std::vector<history_material_base<_T>*> _history_materials;
};
#endif // COMPOSITE_MATERIAL_HISTORY_BASE_BONES_H
