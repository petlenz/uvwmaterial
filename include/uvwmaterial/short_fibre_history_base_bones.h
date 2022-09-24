/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SHORT_FIBRE_HISTORY_BASE_BONES_H
#define SHORT_FIBRE_HISTORY_BASE_BONES_H

#include <cstdlib>

template <typename _T, std::size_t _Dim, typename _Container>
class short_fibre_history_base :
        public history_material_base<_T>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    short_fibre_history_base() = delete;

    short_fibre_history_base(short_fibre_composite_base<_T, _Dim, _Container>* __material);

    virtual ~short_fibre_history_base(){}

    constexpr inline void init(){
        size_type number_of_his_var{0};

        for(auto material : _material->composite_materials()){
            auto history_material = dynamic_cast<history_material_base<_T>*>(material);
            if(history_material){
                _history_materials.push_back(history_material);
                number_of_his_var += history_material->size();
            }
        }
        this->_history.resize(number_of_his_var);
    }

//protected:
    constexpr inline auto set_local_history(){
        size_type iter{0};
        for(auto material : _history_materials){
            const auto his_number{material->size()};
            material->set_history(&this->_history[iter], &this->_history[iter+his_number]);
            iter += his_number;
        }
    }

    constexpr inline auto get_local_history(){
        std::vector<value_type> his;
        his.reserve(this->_history.size());
        for(auto material : _history_materials){
            for(const auto& h : material->history()){
                his.push_back(h);
            }
        }
        this->_history = his;
    }

private:
    short_fibre_composite_base<_T, _Dim, _Container>* _material;
    std::vector<history_material_base<_T>*> _history_materials;
};

template <typename _T, std::size_t _Dim, typename _Container>
short_fibre_history_base<_T, _Dim, _Container>::short_fibre_history_base(short_fibre_composite_base<_T, _Dim, _Container>* __material):
    _material(__material)
{}

#endif // SHORT_FIBRE_HISTORY_BASE_BONES_H
