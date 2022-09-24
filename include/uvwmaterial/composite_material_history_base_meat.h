/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef COMPOSITE_MATERIAL_HISTORY_BASE_MEAT_H
#define COMPOSITE_MATERIAL_HISTORY_BASE_MEAT_H

template <typename _T, std::size_t _Dim, typename _Container>
composite_history_base<_T, _Dim, _Container>::composite_history_base(composite_material_base<_T, _Dim, _Container> * __material):
    _material(__material),
    _history_materials()
{}

template <typename _T, std::size_t _Dim, typename _Container>
composite_history_base<_T, _Dim, _Container>::~composite_history_base(){}

template <typename _T, std::size_t _Dim, typename _Container>
constexpr inline void composite_history_base<_T, _Dim, _Container>::init(){

    for(size_type i{0}; i<_material->number_of_materials(); ++i){
        auto his_material{dynamic_cast<history_material_base<value_type>*>(_material->material(i))};
        if(his_material){
            _history_materials.push_back(his_material);
        }
    }

    auto number_of_his_var{0};
    for(auto his_material : _history_materials){
        number_of_his_var += his_material->size();
    }
    this->_history.resize(number_of_his_var);
}

template <typename _T, std::size_t _Dim, typename _Container>
inline void composite_history_base<_T, _Dim, _Container>::set_local_history(){
    size_type iter{0};
    for(auto his_material : _history_materials){
        const auto his_number{his_material->size()};
        his_material->set_history(&this->_history[iter], &this->_history[iter+his_number]);
        iter += his_number;
    }
}

template <typename _T, std::size_t _Dim, typename _Container>
inline void composite_history_base<_T, _Dim, _Container>::get_local_history(){
    //std::vector<value_type> his;
    //his.reserve(this->_history.size());
    for(auto his_material : _history_materials){
        auto his_com = dynamic_cast<composite_history_base<_T, _Dim, _Container>*>(his_material);
        if(his_com){
            his_com->get_local_history();
        }
    }
    size_type iter{0};
    for(auto his_material : _history_materials){
        for(const auto& h : his_material->history()){
            this->_history[iter++] = h;
            //his.push_back(h);
        }
    }
    //this->_history = his;
}


#endif // COMPOSITE_MATERIAL_HISTORY_BASE_MEAT_H
