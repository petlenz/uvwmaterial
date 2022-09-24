/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef NONLOCAL_COMPOSITE_BASE_BONES_H
#define NONLOCAL_COMPOSITE_BASE_BONES_H

template <typename T, std::size_t Dim, typename Container>
class composite_nonlocal_material_base :
        public nonlocal_material_base<T, Dim>
{
public:
    using value_type = T;
    using size_type = std::size_t;
    using tensor2 = tmech::tensor<value_type, Dim, 2>;
    using nonlocal_key = typename nonlocal_material_base<T, Dim>::key;

    composite_nonlocal_material_base() = delete ;

    composite_nonlocal_material_base(material_base<value_type, Dim, Container> * const __material):
        nonlocal_material_base<T, Dim>(),
        _material(__material)
    {}

    virtual ~composite_nonlocal_material_base(){}

    constexpr inline void init(){
        const auto material_composite{dynamic_cast<composite_material_solid_base<T, Dim, Container>*>(_material)};
        const auto number_of_materials{material_composite->number_of_materials()};

        //set nonlocal domains
        for(size_type i{0}; i<number_of_materials; ++i){
            const auto nonlocal_material{dynamic_cast<nonlocal_material_base<T, Dim>*>(material_composite->material(i))};
            if(nonlocal_material){
                _nonlocal_materials.push_back(nonlocal_material);
            }
        }

        std::map<nonlocal_key, size_type> nonlocal_domains;
        for(auto nonlocal_material : _nonlocal_materials){
            for(const auto& domain : nonlocal_material->nonlocal_domains()){
                nonlocal_domains[domain] += nonlocal_material->number_of_variables(domain);
            }
        }
        for(const auto& [domain, number] : nonlocal_domains){
            this->set_nonlocal_domain(domain, number);
        }

        for(auto nonlocal_material : _nonlocal_materials){
            for(const auto& domain : nonlocal_material->nonlocal_domains()){
                this->_internal_lengths[domain] =  nonlocal_material->internal_length(domain);
                this->_local_transformation[domain] = nonlocal_material->local_transformation(domain);
            }
        }

        _material_composite = dynamic_cast<composite_material_solid_base<T, Dim, Container>*>(_material);
        _composite_history = dynamic_cast<composite_history_base<T, Dim, Container>*>(_material);
    }


    constexpr inline auto update(){
        const auto material_composite{dynamic_cast<composite_material_solid_base<T, Dim, Container>const*>(_material)};
        const auto number_of_materials{material_composite->number_of_materials()};
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

        std::map<nonlocal_key, size_type> iter;
        for(size_type i{0}; i<number_of_materials; ++i){
            const auto nonlocal_material{dynamic_cast<nonlocal_material_base<T,Dim>const*>(material_composite->material(i))};
            if(nonlocal_material){
                for(const auto& domain : nonlocal_material->nonlocal_domains()){
                    if(nonlocal_material->element_marker(domain)){
                        for(size_type j{0}; j<nonlocal_material->receiver_tensor(domain).size(); ++j){
                            this->receiver_tensor(domain)[iter[domain]] = nonlocal_material->receiver_tensor(domain)[j]*material_composite->volume_fraction(i);
                            this->source_tensor(domain)[iter[domain]++] = tmech::dcontract(nonlocal_material->source_tensor(domain)[j], material_composite->strain_concentration_tensors()[i]);
                        }
                        this->element_marker(domain) = true;
                    }
                }
            }
        }
    }


    inline void update_nonlocal_variables(){

        //Transfer strain to local parts
        _material_composite->update_strain();

        //Transfer history to local parts
        if(_composite_history){
            _composite_history->set_local_history();
        }

        //update nonlocal variables
        for(auto nonlocal_material : _nonlocal_materials){
            nonlocal_material->update_nonlocal_variables();
        }

        reset_iter();

        for(auto nonlocal_material : _nonlocal_materials){
            for(const auto& domain : nonlocal_material->nonlocal_domains()){
                const auto& var_res{nonlocal_material->nonlocal_variables(domain)};
                auto& var{this->_nonlocal_variables[domain]};
                auto& iter_local{_iter[domain]};
                for(size_type j{0}; j<var_res.size(); ++j){
                    var[iter_local++] = var_res[j];
                }
            }
        }
//        for(size_type i{0}; i<number_of_materials; ++i){
//            auto nonlocal_material{dynamic_cast<nonlocal_material_base<T, Dim>*>(material_composite->material(i))};
//            if(nonlocal_material){
//                nonlocal_material->update_nonlocal_variables();
//            }
//        }

        //copy nonlocal variables
//        std::unordered_map<nonlocal_key, size_type> iter;
//        for(size_type i{0}; i<number_of_materials; ++i){
//            auto nonlocal_material{dynamic_cast<nonlocal_material_base<T, Dim>*>(material_composite->material(i))};
//            if(nonlocal_material){
//                for(const auto& domain : nonlocal_material->nonlocal_domains()){
//                    const auto& var_res{nonlocal_material->nonlocal_variables(domain)};
//                    auto& var{this->_nonlocal_variables[domain]};
//                    auto& iter_local{iter[domain]};
//                    for(size_type j{0}; j<var_res.size(); ++j){
//                        var[iter_local++] = var_res[j];
//                    }
//                }
//            }
//        }

//        std::unordered_map<std::string, size_type> iter;
//        for(size_type i{0}; i<number_of_materials; ++i){
//            auto nonlocal_material{dynamic_cast<nonlocal_material_base<T,Dim>*>(material_composite->material(i))};
//            if(nonlocal_material){
//                nonlocal_material->update_nonlocal_variables();
//                for(const auto& domain : nonlocal_material->nonlocal_domains()){
//                    for(size_type j{0}; j<nonlocal_material->number_of_variables(domain); ++j){
//                        this->_nonlocal_variables[domain][iter[domain]++] = nonlocal_material->nonlocal_variables(domain)[j];
//                    }
//                }
//            }
//        }
    }

    constexpr inline auto const& nonlocal_materials()const{
        return _nonlocal_materials;
    }

protected:
    constexpr inline auto reset_iter(){
        for(auto & [domain, iter] : _iter){
            iter = 0;
        }
    }

    inline void set_internal_nonlocal_variables(){
        //const auto material_composite{dynamic_cast<composite_material_base<T, Dim, Container>*>(_material)};
        //const auto number_of_materials{material_composite->number_of_materials()};

        std::map<nonlocal_key, size_type> iter;
        for(auto nonlocal_material : _nonlocal_materials){
            for(const auto& domain : nonlocal_material->nonlocal_domains()){
                const auto number_nonlocal_variables{nonlocal_material->number_of_variables(domain)};
                std::vector<value_type> temp(number_nonlocal_variables);
                for(size_type j{0}; j<number_nonlocal_variables; ++j){
                    temp[j] = this->nonlocal_variables(domain)[iter[domain]++];
                }
                nonlocal_material->nonlocal_variables(domain) = temp;
            }
        }

//        for(size_type i{0}; i<number_of_materials; ++i){
//            auto nonlocal_material{dynamic_cast<nonlocal_material_base<T,Dim>*>(material_composite->material(i))};
//            if(nonlocal_material){
//                for(const auto& domain : nonlocal_material->nonlocal_domains()){
//                    const auto number_nonlocal_variables{nonlocal_material->number_of_variables(domain)};
//                    std::vector<value_type> temp(number_nonlocal_variables);
//                    for(size_type j{0}; j<number_nonlocal_variables; ++j){
//                        temp[j] = this->nonlocal_variables(domain)[iter[domain]++];
//                    }
//                    nonlocal_material->nonlocal_variables(domain) = temp;
//                }
//            }
//        }
    }

private:
    std::unordered_map<nonlocal_key, size_type> _iter;
    material_base<T, Dim, Container> * const _material;
    std::vector<nonlocal_material_base<T, Dim>*> _nonlocal_materials;
    composite_history_base<T, Dim, Container>* _composite_history;
    composite_material_solid_base<T, Dim, Container>* _material_composite;
};

#endif // NONLOCAL_COMPOSITE_BASE_BONES_H
