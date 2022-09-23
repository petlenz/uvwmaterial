#ifndef SHORT_FIBRE_NONLOCAL_BASE_BONES_H
#define SHORT_FIBRE_NONLOCAL_BASE_BONES_H


#include "small_strain_material_base_bones.h"

template <typename _T, std::size_t _Dim, typename _Container>
class short_fibre_nonlocal_base :
        public nonlocal_material_base<_T, _Dim>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;
    using tensor2 = tmech::tensor<value_type, _Dim, 2>;
    using nonlocal_key = typename nonlocal_material_base<_T, _Dim>::key;

    short_fibre_nonlocal_base() = delete;

    short_fibre_nonlocal_base(short_fibre_composite_base<_T, _Dim, _Container>* __material):
        nonlocal_material_base<_T, _Dim>(),
        _material(__material),
        _nonlocal_materials(),
        _all_nonlocal_materials(),
        _all_materials(),
        _composite_materials(),
        _iter(),
        _test()
    {}

    virtual ~short_fibre_nonlocal_base(){}

    constexpr inline void init(){
        for(auto material : _material->composite_materials()){
            const auto composite_nonlocal{dynamic_cast<composite_nonlocal_material_base<_T, _Dim, _Container>*>(material)};
            if(composite_nonlocal){
                for(auto mat : composite_nonlocal->nonlocal_materials()){
                    _all_nonlocal_materials.push_back(mat);
                }
            }
            const auto nonlocal_material{dynamic_cast<nonlocal_material_base<_T, _Dim>*>(material)};
            if(nonlocal_material){
                _nonlocal_materials.push_back(nonlocal_material);
            }
        }


        for(auto material : _material->composite_materials()){
            _composite_materials.push_back(dynamic_cast<material_base<_T, _Dim, _Container>*>(material));
            for(auto mat : material->materials()){
                _all_materials.push_back(mat);
            }
        }

        //set nonlocal domains
        std::map<nonlocal_key, size_type> nonlocal_domains;
        for(const auto* mat : _nonlocal_materials){
            for(const auto& domain : mat->nonlocal_domains()){
                nonlocal_domains[domain] += mat->number_of_variables(domain);
            }
        }
        for(const auto& [domain, number] : nonlocal_domains){
            this->set_nonlocal_domain(domain, number);
        }
        for(const auto* mat : _nonlocal_materials){
            for(const auto& domain : mat->nonlocal_domains()){
                this->_internal_lengths[domain] =  mat->internal_length(domain);
                this->_local_transformation[domain] = mat->local_transformation(domain);
            }
        }

        //iter
        for(const auto&  nonlocal_domain : this->_nonlocal_domains){
            _iter[nonlocal_domain] = 0;
        }

        //test _nonlocal_materials
        for(const auto* mat : _all_nonlocal_materials){
            for(const auto& [domain, var] : mat->nonlocal_variables()){
                _test.emplace_back(std::make_tuple(&_iter[domain], &(this->_nonlocal_variables.find(domain)->second), &var));
            }
        }
    }

    constexpr inline auto update(){
        for(auto& [_, vec_tensor] : this->_source){
            for(auto& tensor : vec_tensor)
                tensor.fill(0);
        }

        for(auto& [_, vec_tensor] : this->_receiver){
            for(auto& tensor : vec_tensor)
                tensor.fill(0);
        }

        for(auto& [_, marker] : this->_element_marker){
            marker = false;
        }


        reset_iter();

        const auto& weights{_material->weights()};
        for(size_type i{0}; i<_nonlocal_materials.size(); ++i){
            const auto* nonlocal_material{_nonlocal_materials[i]};
            for(const auto& domain : nonlocal_material->nonlocal_domains()){
                if(nonlocal_material->element_marker(domain)){
                    for(size_type j{0}; j<nonlocal_material->receiver_tensor(domain).size(); ++j){
                        this->receiver_tensor(domain)[_iter[domain]] = weights[i]*nonlocal_material->receiver_tensor(domain)[j];
                        this->source_tensor(domain)[_iter[domain]++] = nonlocal_material->source_tensor(domain)[j];
                    }
                    this->element_marker(domain) = true;
                }
            }
        }
    }


    inline void update_nonlocal_variables(){
        const auto composite_history{dynamic_cast<short_fibre_history_base<_T, _Dim, _Container>*>(_material)};

        //Transfer history to local parts
        if(composite_history){
            composite_history->set_local_history();
        }

        //Transfer strain to local parts
        _material->update_strain();

        //update nonlocal variabels _nonlocal_materials
        for(auto nonlocal_material : _nonlocal_materials){
            nonlocal_material->update_nonlocal_variables();
        }

        reset_iter();

        //<ptr_iter, ptr_to_variable_rec, ptr_to_variable_src>
        for(auto [iter, ptr_rec, ptr_src] : _test){
            const auto& src{*ptr_src};
            auto& rec{*ptr_rec};
            for(size_type j{0}; j<src.size(); ++j){
                rec[*iter++] = src[j];
            }
        }

//        for(auto nonlocal_material : _nonlocal_materials){
//            for(const auto& [domain, var_res] : nonlocal_material->nonlocal_variables()){
//                auto& var{this->_nonlocal_variables[domain]};
//                auto& iter_local{_iter[domain]};
//                for(size_type j{0}; j<var_res.size(); ++j){
//                    var[iter_local++] = var_res[j];
//                }
//            }
//        }
    }

protected:

    constexpr inline auto reset_iter(){
        for(auto & [domain, iter] : _iter){
            iter = 0;
        }
    }

    inline void set_internal_nonlocal_variables(){
        reset_iter();
        //size_type iter{0};
        //std::map<nonlocal_key, size_type> iter;
        for(auto material : _nonlocal_materials){
            for(const auto& domain : material->nonlocal_domains()){
                const auto number_nonlocal_variables{material->number_of_variables(domain)};
                std::vector<value_type> temp(number_nonlocal_variables);
                for(size_type j{0}; j<number_nonlocal_variables; ++j){
                    temp[j] = this->nonlocal_variables(domain)[_iter[domain]++];
                }
                material->nonlocal_variables(domain) = temp;
            }
        }
    }

private:
    short_fibre_composite_base<_T, _Dim, _Container>* _material;
    std::vector<nonlocal_material_base<_T,_Dim>*> _nonlocal_materials;
    std::vector<nonlocal_material_base<_T, _Dim>*> _all_nonlocal_materials;
    std::vector<material_base<_T, _Dim, _Container>*> _all_materials;
    std::vector<material_base<_T, _Dim, _Container>*> _composite_materials;
    std::unordered_map<nonlocal_key, size_type> _iter;
    std::vector<std::tuple<size_type*, std::vector<value_type>*, std::vector<value_type> const*>> _test;
    //std::map<std::string, std::string> _matrix_and_fibre_nonlocal_domain;
};

#endif // SHORT_FIBRE_NONLOCAL_BASE_BONES_H
