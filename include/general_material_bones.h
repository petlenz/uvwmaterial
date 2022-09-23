#ifndef GENERAL_MATERIAL_BONES_H
#define GENERAL_MATERIAL_BONES_H

#include <vector>
#include <memory>
#include <math.h>


template <typename T, std::size_t Dim, typename Container>
class general_material
{
    template<typename Type>
    using vector = std::vector<std::unique_ptr<Type>>;

public:
    using size_type = std::size_t;
    using value_type = T;
    using tensor2map = tmech::adaptor<T, Dim, 2, tmech::full<Dim>>;
    using tensor1map = tmech::adaptor<T, Dim, 1, tmech::full<Dim>>;

    general_material():
        material(),
        base_materials(),
        propagation_laws(),
        yield_functions(),
        state_functions(),
        inclusions(),
        eshelby_tensors(),
        r_o_m_kernals(),
        mean_field_composite_kernals(),
        curing_functions()
    {}

    template<typename Json>
    constexpr inline auto make_material(Json const& json){
        material = make_material_detail(json);
    }

    constexpr inline auto clear(){
        material.clear();
        base_materials.clear();
        propagation_laws.clear();
        yield_functions.clear();
        state_functions.clear();
        inclusions.clear();
        eshelby_tensors.clear();
    }

    constexpr inline auto* get(){
        return material.get();
    }

    constexpr inline auto* get()const{
        return material.get();
    }

private:
    template<typename Json>
    constexpr inline std::unique_ptr<material_base<T, Dim, Container>> make_material_detail(Json const& json){
        const std::string name{json["name"]};
        if(json["name"] == "linear_elasticity"){
            return make_linear_elasticity<T, Dim, Container>(json["parameter"]["E"], json["parameter"]["nu"]);
        }/*else if(json["name"] == "linear_thermo_elasticity"){
            return make_linear_thermo_elasticity<T, Dim, Container>(json["parameter"]["E"], json["parameter"]["nu"], json["parameter"]["alpha"]);
        }*/else if(json["name"] == "small_strain_isotropic_damage"){
            //base material
            base_materials.emplace_back(make_material_detail(json["base_material"]));

            //state function
            state_functions.emplace_back(make_state_function(json["state_function"]));
            state_functions.back().get()->set_base_material(make_solid_material(base_materials.back().get()));

            //yield function
            yield_functions.emplace_back(make_yield_function(json["yield_function"]));
            yield_functions.back().get()->base_material(make_solid_material(base_materials.back().get()));
            yield_functions.back().get()->state_function(0,state_functions.back().get());

            //propagation law
            propagation_laws.emplace_back(make_propagation_law(json["propagation_law"]));

            //final material
            auto mat_base = make_small_strain_material_isotropic_damage<T, Dim, Container>();
            auto mat = make_small_strain_material_isotropic_damage(mat_base.get());
            mat->base_material(make_small_strain_material(base_materials.back().get()));
            mat->yield_function(make_yield_function_isotropic_damage_strain_based_base(yield_functions.back().get()));
            mat->propagation_law(propagation_laws.back().get());
            return mat_base;
        }
        /*else if(json["name"] == "small_strain_plasticity_single_function"){
            //base material
            base_materials.emplace_back(make_material_detail(json["base_material"]));

            //hardening stress
            propagation_laws.emplace_back(make_propagation_law(json["hardening_stress"]));

            //yield function
            yield_functions.emplace_back(make_yield_function(json["yield_function"]));
            yield_functions.back().get()->set_base_material(base_materials.back().get());
            auto stress_yield{make_yield_single_function_stress_base(yield_functions.back().get())};
            stress_yield->set_hardening_stress(propagation_laws.back().get());

            //final material
            auto mat = make_small_strain_material_plasticity_single_yield_function<T, Dim, Container>();
            mat->set_base_material(make_small_strain_material(base_materials.back().get()));
            mat->set_yield_function(stress_yield);
            mat->set_parameter(json["parameter"]["G"]);
            return mat;
        }*/
        if(json["name"] == "small_strain_rule_of_mixture_composite"){
            //final material
            auto mat = make_small_strain_rule_of_mixture<T, Dim, Container>();
            r_o_m_kernals.push_back(make_rule_of_mixture_kernal(json));
            mat->set_kernal(r_o_m_kernals.back().get());
            for(const auto material : json["materials"]){
                base_materials.push_back(make_material_detail(material));
                mat->push_back(material["volume_fraction"], make_small_strain_material(base_materials.back().get()));
            }
            return mat;
        }

        if(json["name"] == "small_strain_rule_of_mixture_history"){
            //final material
            auto mat = make_small_strain_rule_of_mixture_history<T, Dim, Container>();
            r_o_m_kernals.push_back(make_rule_of_mixture_kernal(json));
            mat->set_kernal(r_o_m_kernals.back().get());
            for(const auto material : json["materials"]){
                base_materials.push_back(make_material_detail(material));
                mat->push_back(material["volume_fraction"], make_small_strain_material(base_materials.back().get()));
            }
            return mat;
        }

        if(json["name"] == "small_strain_mean_field_composite"){
            //final material
            auto mat = make_small_strain_mean_field_composite<T, Dim, Container>();
            mean_field_composite_kernals.push_back(make_mean_field_composite_kernal(json));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            return mat;
        }

        if(json["name"] == "small_strain_mean_field_composite_history"){
            //final material
            auto mat = make_small_strain_mean_field_composite_history<T, Dim, Container>();
            mean_field_composite_kernals.push_back(make_mean_field_composite_kernal(json));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            return mat;
        }

        if(json["name"] == "small_strain_mean_field_composite_damage_nonlocal"){
            //final material
            auto mat = make_small_strain_mean_field_composite_damage_nonlocal<T, Dim, Container>();
            mean_field_composite_kernals.push_back(make_mean_field_composite_kernal(json));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            return mat;
        }

        if(json["name"] == "small_strain_mean_field_composite_plastic"){
            //final material
            auto mat = make_small_strain_mean_field_composite_history<T, Dim, Container>();
            mean_field_composite_kernals.push_back(make_mean_field_composite_kernal(json));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            //kernal nonlinear_mean_field
            //"type":"gradient"
            //"type":"newton"
            return mat;
        }

        if(json["name"] == "small_strain_isotropic_nonlocal_damage"){
            //equivalent_strain
            //damage_variable

            //base material
            base_materials.emplace_back(make_material_detail(json["base_material"]));

            //state function
            state_functions.emplace_back(make_state_function(json["state_function"]));
            state_functions.back().get()->set_base_material(static_cast<solid_material_base<T,Dim, Container>*>(base_materials.back().get()));

            //yield function
            yield_functions.emplace_back(make_yield_function(json["yield_function"]));
            yield_functions.back().get()->base_material(make_solid_material(base_materials.back().get()));
            yield_functions.back().get()->state_function(0, state_functions.back().get());

            //propagation law
            propagation_laws.emplace_back(make_propagation_law(json["propagation_law"]));

            //final material
            auto mat_base = make_small_strain_material_isotropic_nonlocal_damage<T, Dim, Container>();
            auto mat = make_small_strain_material_isotropic_nonlocal_damage(mat_base.get());
            mat->set_base_material(make_small_strain_material(base_materials.back().get()));
            mat->set_yield_function(make_yield_function_isotropic_damage_strain_based_base(yield_functions.back().get()));
            mat->set_propagation_law(propagation_laws.back().get());
            mat->average_type() = json["average_type"];
            mat->set_nonlocal_domain(json["nonlocal_domain"], 1);
            mat->internal_length(json["nonlocal_domain"]) = json["internal_length"].template get<std::vector<T>>();
            if(json.contains("local_transformation")){
                std::vector<T> temp = json["local_transformation"].template get<std::vector<T>>();
                mat->local_transformation(json["nonlocal_domain"]) = tensor2map(&temp[0]);
            }
            return mat_base;
        }

        if(json["name"] == "small_strain_short_fibre_composite"){
            auto mat = make_small_strain_short_fibre<T, Dim, Container>();
            std::vector<T> fot = json["fibre_orientation_tensor"].template get<std::vector<T>>();
            mat.get()->set_fibre_orientation_distribution_tensor(tensor2map(&fot[0]));
            //init matrix B and integration points
            auto short_fibre_mat{dynamic_cast<short_fibre_composite_base<T, Dim, Container>*>(mat.get())};
            short_fibre_mat->init();
            const auto& directions{short_fibre_mat->directions()};
            //loop over integration points
            for(size_type i{0}; i<mat.get()->number_of_composites(); ++i){
                base_materials.emplace_back(make_material_detail(json["composite"]));
                if(json["composite"]["inclusions"]["inclusion1"]["eshelby_tensor"]["type"] == "ellipsoid"){
                    dynamic_cast<eshelby_tensor_solid_cylindrical_ellipsoid<T, Dim>*>(eshelby_tensors.back().get())->direction() = directions[i];
                }else{
                    throw std::runtime_error("general_material::make_material_detail(): no matching eshelby tensor");
                }
                mat.get()->push_back(dynamic_cast<composite_material_base<T, Dim, Container>*>(base_materials.back().get()));
            }
            return mat;
        }

        if(json["name"] == "small_strain_short_fibre_history_composite"){
            auto mat = make_small_strain_short_fibre_history<T, Dim, Container>();
            std::vector<T> fot = json["fibre_orientation_tensor"].template get<std::vector<T>>();
            tensor2map fot_t(&fot[0]);
            std::cout<<fot_t<<std::endl;
            mat.get()->set_fibre_orientation_distribution_tensor(fot_t);
            mat.get()->integration_precision() = json["integration_precision"];
            //init matrix B and integration points
            auto short_fibre_mat{dynamic_cast<short_fibre_composite_base<T, Dim, Container>*>(mat.get())};
            short_fibre_mat->init();
            const auto& directions{short_fibre_mat->directions()};
            //loop over integration points
            for(size_type i{0}; i<mat.get()->number_of_composites(); ++i){
                //make composite material
                base_materials.emplace_back(make_material_detail(json["composite"]));


                if(json["composite"]["inclusions"]["inclusion1"]["eshelby_tensor"]["type"] == "ellipsoid"){
                    dynamic_cast<eshelby_tensor_solid_cylindrical_ellipsoid<T, Dim>*>(eshelby_tensors.back().get())->direction() = directions[i];
                }else{
                    throw std::runtime_error("general_material::make_material_detail(): no matching eshelby tensor");
                }
                if(json["composite"]["inclusions"]["inclusion1"]["material"].contains("state_function")){
                    if constexpr (Dim == 3){
                        if(json["composite"]["inclusions"]["inclusion1"]["material"]["state_function"]["name"] == "vector_strain_state_function"){
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*2].get())->set_direction_vector(directions[i]);
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*2].get())->set_tangential_vector(directions[i]);
                            std::cout<<directions[i]<<std::endl;
                        }
                    }else{
                        throw std::runtime_error("general_material::short_fibre_history_composite()");
                    }
                }
                mat.get()->push_back(dynamic_cast<composite_material_base<T, Dim, Container>*>(base_materials.back().get()));
            }
            return mat;
        }

        if(json["name"] == "small_strain_short_fibre_damage_nonlocal"){
            auto mat = make_small_strain_short_fibre_damage_nonlocal<T, Dim, Container>();
            std::vector<T> fot = json["fibre_orientation_tensor"].template get<std::vector<T>>();
            mat.get()->set_fibre_orientation_distribution_tensor(tensor2map(&fot[0]));
            mat.get()->integration_precision() = json["integration_precision"];
            //init matrix B and integration points
            auto short_fibre_mat{dynamic_cast<short_fibre_composite_base<T, Dim, Container>*>(mat.get())};
            short_fibre_mat->init();
            const auto& directions{short_fibre_mat->directions()};
            //loop over integration points
            for(size_type i{0}; i<mat.get()->number_of_composites(); ++i){
                //local copy to change nonlocal domain
                Json composite = json["composite"];

                //set inclusion and matrix nonlocal domain
                if(json["composite"]["inclusions"]["inclusion1"]["material"].contains("nonlocal_domain")){
                    const size_type nonlocal = json["composite"]["inclusions"]["inclusion1"]["material"]["nonlocal_domain"];
                    composite["inclusions"]["inclusion1"]["material"]["nonlocal_domain"] = std::stoi(std::to_string(nonlocal) + std::to_string(i));
                }
                if(json["composite"]["matrix"]["material"].contains("nonlocal_domain")){
                    const size_type nonlocal = json["composite"]["matrix"]["material"]["nonlocal_domain"];
                    composite["matrix"]["material"]["nonlocal_domain"] = std::stoi(std::to_string(nonlocal) + std::to_string(i));
                }

                base_materials.emplace_back(make_material_detail(composite));
                //set direction of eshelby's tensor
                if(json["composite"]["inclusions"]["inclusion1"]["eshelby_tensor"]["type"] == "ellipsoid"){
                    dynamic_cast<eshelby_tensor_solid_cylindrical_ellipsoid<T, Dim>*>(eshelby_tensors.back().get())->direction() = directions[i];
                }else{
                    throw std::runtime_error("general_material::make_material_detail(): no matching eshelby tensor");
                }
                //set directions of state function
                if(json["composite"]["inclusions"]["inclusion1"]["material"].contains("state_function")){
                    if constexpr (Dim == 3){
                        if(json["composite"]["inclusions"]["inclusion1"]["material"]["state_function"]["name"] == "vector_strain_state_function"){
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*2].get())->set_direction_vector(directions[i]);
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*2].get())->set_tangential_vector(directions[i]);
                        }
                    }else{
                        throw std::runtime_error("general_material::short_fibre_history_composite()");
                    }
                }

                //set local transformation tensor in terms of the direction vectors
                auto composit{dynamic_cast<composite_material_base<T, Dim, Container>*>(base_materials.back().get())};
                for(size_type ii{0}; ii<composit->number_of_materials(); ++ii){
                    auto nonlocal_mat = dynamic_cast<nonlocal_material_base<T, Dim>*>(composit->material(ii));
                    if(nonlocal_mat){
                        for(auto const& domain : nonlocal_mat->nonlocal_domains()){
                            if constexpr (Dim == 3){
                                nonlocal_mat->local_transformation(domain) = rotation_direction(directions[i]);
                                //std::cout<<tmech::otimes(directions[i], directions[i])<<std::endl;
                                //std::cout<<rotation_direction(tmech::abs(directions[i]))<<std::endl;
                                //std::cout<<rotation_euler(std::acos(directions[i](2)), std::acos(directions[i](1)), std::acos(directions[i](0)))<<std::endl;
                            }else{
                                throw std::runtime_error("general_material::short_fibre_history_composite()");
                            }
                        }
                    }
                }

                mat.get()->push_back(composit);
            }
            return mat;
        }

        if(json["name"] == "small_strain_short_fibre_interface_damage_nonlocal"){
            auto mat = make_small_strain_short_fibre_damage_nonlocal<T, Dim, Container>();
            std::vector<T> fot = json["fibre_orientation_tensor"].template get<std::vector<T>>();
            mat.get()->set_fibre_orientation_distribution_tensor(tensor2map(&fot[0]));
            mat.get()->integration_precision() = json["integration_precision"];
            //init matrix B and integration points
            auto short_fibre_mat{dynamic_cast<short_fibre_composite_base<T, Dim, Container>*>(mat.get())};
            short_fibre_mat->init();
            const auto& directions{short_fibre_mat->directions()};
            //loop over integration points
            for(size_type i{0}; i<mat.get()->number_of_composites(); ++i){
                //local copy to change nonlocal domain
                Json composite = json["composite"];

                //set inclusion and matrix nonlocal domain
                if(json["composite"]["inclusions"]["inclusion1"]["material"].contains("nonlocal_domain")){
                    auto& material{composite["composite"]["inclusions"]["inclusion1"]["material"]};
                    //fibre
                    if(material["inclusions"]["inclusion1"]["material"].contains("nonlocal_domain")){
                        const size_type nonlocal = material["inclusions"]["inclusion1"]["material"]["nonlocal_domain"];
                        material["inclusions"]["inclusion1"]["material"]["nonlocal_domain"] = std::stoi(std::to_string(nonlocal) + std::to_string(i));
                    }
                    //interface
                    if(material["matrix"]["material"].contains("nonlocal_domain")){
                        const size_type nonlocal = material["material"]["material"]["nonlocal_domain"];
                        material["matrix"]["material"]["nonlocal_domain"] = std::stoi(std::to_string(nonlocal) + std::to_string(i));
                    }
                }

                if(json["composite"]["matrix"]["material"].contains("nonlocal_domain")){
                    const size_type nonlocal = json["composite"]["matrix"]["material"]["nonlocal_domain"];
                    composite["matrix"]["material"]["nonlocal_domain"] = std::stoi(std::to_string(nonlocal) + std::to_string(i));
                }


                base_materials.emplace_back(make_material_detail(composite));
                //set direction of eshelby's tensor
                if(json["composite"]["inclusions"]["inclusion1"]["eshelby_tensor"]["type"] == "ellipsoid"){
                    dynamic_cast<eshelby_tensor_solid_cylindrical_ellipsoid<T, Dim>*>(eshelby_tensors.back().get())->direction() = directions[i];
                }else{
                    throw std::runtime_error("general_material::make_material_detail(): no matching eshelby tensor");
                }
                //set directions of state function in fibre
                if(json["composite"]["inclusions"]["inclusion1"]["material"]["inclusions"]["inclusion1"]["material"].contains("state_function")){
                    if constexpr (Dim == 3){
                        if(json["composite"]["inclusions"]["inclusion1"]["material"]["inclusions"]["inclusion1"]["material"]["state_function"]["name"] == "vector_strain_state_function"){
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*3].get())->set_direction_vector(directions[i]);
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*3].get())->set_tangential_vector(directions[i]);
                        }
                    }else{
                        throw std::runtime_error("general_material::short_fibre_history_composite()");
                    }
                }

                //set directions of state function in fibre
                if(json["composite"]["inclusions"]["inclusion1"]["material"]["matrix"]["material"].contains("state_function")){
                    if constexpr (Dim == 3){
                        if(json["composite"]["inclusions"]["inclusion1"]["material"]["matrix"]["material"]["state_function"]["name"] == "vector_strain_state_function"){
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*3+1].get())->set_direction_vector(directions[i]);
                            dynamic_cast<vector_strain_state_function<T, Dim, Container>*>(state_functions[i*3+1].get())->set_tangential_vector(directions[i]);
                        }
                    }else{
                        throw std::runtime_error("general_material::short_fibre_history_composite()");
                    }
                }

                //set local transformation tensor in terms of the direction vectors
                auto composit{dynamic_cast<composite_material_base<T, Dim, Container>*>(base_materials.back().get())};
                for(size_type ii{0}; ii<composit->number_of_materials(); ++ii){
                    auto nonlocal_mat = dynamic_cast<nonlocal_material_base<T, Dim>*>(composit->material(ii));
                    if(nonlocal_mat){
                        for(auto const& domain : nonlocal_mat->nonlocal_domains()){
                            if constexpr (Dim == 3){
                                nonlocal_mat->local_transformation(domain) = rotation_direction(directions[i]);
                                //std::cout<<tmech::otimes(directions[i], directions[i])<<std::endl;
                                //std::cout<<rotation_direction(tmech::abs(directions[i]))<<std::endl;
                                //std::cout<<rotation_euler(std::acos(directions[i](0))*180.0/M_PI, std::acos(directions[i](1))*180.0/M_PI, std::acos(directions[i](2))*180.0/M_PI)<<std::endl;
                            }else{
                                throw std::runtime_error("general_material::short_fibre_history_composite()");
                            }
                        }
                    }
                }

                mat.get()->push_back(composit);
            }
            return mat;
        }


        if(json["name"] == "small_strain_plasticity_single_function"){
            //base material
            base_materials.emplace_back(make_material_detail(json["base_material"]));

            //hardening stress
            propagation_laws.emplace_back(make_propagation_law(json["hardening_stress"]));

            //yield function
            yield_functions.emplace_back(make_yield_function(json["yield_function"]));
            yield_functions.back().get()->base_material(make_solid_material(base_materials.back().get()));
            auto stress_yield{make_yield_function_stress_base(yield_functions.back().get())};
            stress_yield->set_hardening_stress(propagation_laws.back().get());

            //final material
            auto mat = make_small_strain_material_plasticity_single_yield_function<T, Dim, Container>();
            mat->set_base_material(make_small_strain_material(base_materials.back().get()));
            mat->set_yield_function(stress_yield);
            mat->set_parameter(json["parameter"]["G"]);
            return mat;
        }

        if(json["name"] == "linear_elastic_incremental"){
            auto mat = std::make_unique<incremental_linear_elasticity<T, Dim, Container>>();
            mat->set_parameter(json["parameter"]["E"]);
            mat->set_parameter(json["parameter"]["nu"]);
            return mat;
        }


        if(json["name"] == "small_strain_incremental_plasticity_single_function"){
            //base material
            base_materials.emplace_back(make_material_detail(json["base_material"]));

            //hardening stress
            propagation_laws.emplace_back(make_propagation_law(json["hardening_stress"]));

            //yield function
            yield_functions.emplace_back(make_yield_function(json["yield_function"]));
            yield_functions.back().get()->base_material(make_solid_material(base_materials.back().get()));
            auto stress_yield{make_yield_function_stress_base(yield_functions.back().get())};
            stress_yield->set_hardening_stress(propagation_laws.back().get());

            //final material
            auto mat = std::make_unique<small_strain_incremental_plasticity_single_yield_function<T, Dim, Container>>();
            mat->set_base_material(make_small_strain_material(base_materials.back().get()));
            mat->set_yield_function(stress_yield);
            mat->set_parameter(json["parameter"]["G"]);
            return mat;
        }

        //finite strain material
        if(json["name"] == "saint_venant_kirchhoff"){
            auto mat = std::make_unique<saint_venant_kirchhoff<T, Dim, Container>>();
            mat->set_parameter(json["parameter"]["K"]);
            mat->set_parameter(json["parameter"]["G"]);
            return mat;
        }

        if(json["name"] == "incremental_saint_venant_kirchhoff"){
            auto mat = std::make_unique<incremental_saint_venant_kirchhoff<T, Dim, Container>>();
            mat->set_parameter(json["parameter"]["K"]);
            mat->set_parameter(json["parameter"]["G"]);
            return mat;
        }

        if(json["name"] == "hencky_strains_linear_elastic"){
            auto mat = std::make_unique<hencky_strains_linear_elastic<T, Dim, Container>>();
            mat->set_parameter(json["parameter"]["K"]);
            mat->set_parameter(json["parameter"]["G"]);
            return mat;
        }

        if(json["name"] == "hencky_strains_linear_thermo_elastic"){
            auto mat = std::make_unique<hencky_strains_linear_thermo_elastic<T, Dim, Container>>();
            mat->set_parameter(json["parameter"]["K"]);
            mat->set_parameter(json["parameter"]["G"]);
            mat->set_parameter(json["parameter"]["alpha"]);
            return mat;
        }

        if(json["name"] == "hencky_strains_linear_thermo_chemo_elastic"){
            auto mat = std::make_unique<hencky_strains_linear_thermo_chemo_elastic<T, Dim, Container>>();
            mat->set_parameter(json["parameter"]["K"]);
            mat->set_parameter(json["parameter"]["G"]);
            mat->set_parameter(json["parameter"]["alpha"]);
            mat->set_parameter(json["parameter"]["beta"]);
            return mat;
        }

        if(json["name"] == "finite_strain_mean_field_composite"){
            auto mat = std::make_unique<finite_strain_mean_field_composite<T, Dim, Container>>();
            mean_field_composite_kernals.push_back(make_mean_field_composite_kernal(json));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            mat->set_parameter(json["parameter"]["tol"]);
            mat->set_parameter(json["parameter"]["max_iter"]);
            return mat;
        }

        if(json["name"] == "finite_strain_mean_field_composite_explicit"){
            auto mat = std::make_unique<finite_strain_mean_field_composite_explicit<T, Dim, Container>>();
            mean_field_composite_kernals.push_back(make_mean_field_composite_kernal(json));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            //mat->set_parameter(json["parameter"]["volume_fraction"]);
            return mat;
        }

        if(json["name"] == "finite_strain_mean_field_composite_elastic"){
            auto mat = std::make_unique<finite_strain_mean_field_composite_elastic<T, Dim, Container>>();
            mean_field_composite_kernals.push_back(make_mean_field_composite_kernal(json));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            return mat;
        }


        //thermal conductivity
        if(json["name"] == "linear_isotropic_conductivity"){
            auto mat = std::make_unique<linear_isotropic_thermal_conductivity<T, Dim, Container>>();
            mat->set_parameter(json["parameter"]["k"]);
            return mat;
        }

        if(json["name"] == "conductivity_mean_field_composite"){
            //final material
            auto mat = std::make_unique<mean_field_composite_conductivity_material<T, Dim, Container>>();
            auto kernal = make_mean_field_composite_kernal(json);
            kernal->set_material(mat.get());
            mean_field_composite_kernals.push_back(std::move(kernal));
            mat->kernal(mean_field_composite_kernals.back().get());
            make_mean_field_composite_material(mat.get(), json);
            return mat;
        }

        if(json["name"] == "gamm2022"){
            auto mat = std::make_unique<gamm2022<T,Dim,Container>>();
            auto mat_solid = make_material_detail(json["solid_material"]);
            auto mat_thermal = make_material_detail(json["thermal_material"]);
            auto cuing_func = make_curing_function(json["curing_function"]);
            mat->set_solid_material(mat_solid.get());
            mat->set_curing_function(cuing_func.get());
            mat->set_thermo_material(mat_thermal.get());

            base_materials.push_back(std::move(mat_solid));
            base_materials.push_back(std::move(mat_thermal));
            curing_functions.push_back(std::move(cuing_func));

            mat->set_parameter(json["parameter"]["n"]);
            mat->set_parameter(json["parameter"]["cd"]);
            mat->set_parameter(json["parameter"]["rho"]);
            mat->set_parameter(json["parameter"]["Ht"]);
            mat->set_parameter(json["parameter"]["Rn"]);
            mat->set_parameter(json["parameter"]["rho_inc"]);
            mat->set_parameter(json["parameter"]["cd_inc"]);
            return mat;
        }

        throw std::runtime_error("general_material::make_material_detail(): no matching material");
        return nullptr;
    }

    template<typename Json>
    constexpr inline std::unique_ptr<kinetic_model_thermoset_base<T>> make_curing_function(Json const& json){
        if(json["name"] == "autocatalytic_reaction"){
            auto func = std::make_unique<autocatalytic_reaction<T>>();
            auto const& parameter{json["parameter"]};
            func->set_parameter(parameter["A"],parameter["E"],parameter["n"],parameter["m"]);
            return func;
        }
    }

    template<typename Json>
    constexpr inline std::unique_ptr<rule_of_mixture_kernal_base<T, Dim, Container>> make_rule_of_mixture_kernal(Json const& json){
        if(json["kernal"] == "voigt"){
            return make_rule_of_mixture_kernal_voigt<T, Dim, Container>();
        }else if(json["kernal"] == "reuss"){
            return make_rule_of_mixture_kernal_reuss<T, Dim, Container>();
        }else if(json["kernal"] == "vrh"){
            return make_rule_of_mixture_kernal_vrh<T, Dim, Container>();
        }else{
            throw std::runtime_error("general_material::small_strain_rule_of_mixture(): no matching kernal");
        }
    }

    template<typename Json>
    constexpr inline std::unique_ptr<mean_field_composite_kernal_base<T, Dim, Container>> make_mean_field_composite_kernal(Json const& json){
        const std::string name{json["kernal"]};
        if(name == "mori_tanaka"){
            return make_mean_field_mori_tanka_kernal<T, Dim, Container>();
        }

        if(name == "dilute"){
            return make_mean_field_dilute_kernal<T, Dim, Container>();
        }

        if(name == "scs"){
            return make_mean_field_scs_kernal<T, Dim, Container>();
        }

        if(name == "finite_strain_dilute"){
            return std::make_unique<mean_field_composite_solid_finite_strain_dilute_kernal<T, Dim, Container>>();
        }

        if(name == "finite_strain_mori_tanaka"){
            return std::make_unique<mean_field_composite_solid_finite_strain_mori_tanaka_kernal<T, Dim, Container>>();
        }

        if(name == "conductivity_dilute"){
            return std::make_unique<mean_field_composite_conductivity_dilute_kernal<T, Dim, Container>>();
        }

        if(name == "conductivity_mori_tanaka"){
            return std::make_unique<mean_field_composite_conductivity_mori_tanaka_kernal<T, Dim, Container>>();
        }

        throw std::runtime_error("general_material::mean_field_composite_kernal_base(): no matching kernal");
        return nullptr;
    }

    template<typename Json>
    constexpr inline std::unique_ptr<yield_function_base<T, Dim, Container>> make_yield_function(Json const& json){
        if(json["name"] == "strain_based_damage"){
            auto yield_func{make_isotropic_damage_strain_based_yield_function<T, Dim, Container>()};
            yield_func->critical_value() = json["parameter"]["crit_val"];
            return yield_func;
        }else if(json["name"] == "j2_function"){
            auto yield_func{make_j2_yield_function<T, Dim, Container>()};
            yield_func->resize(1);
            yield_func->critical_value(0, json["parameter"]["crit_val"]);
            yield_func->set_parameter(json["parameter"]["G"]);
            state_functions.emplace_back(make_von_mises_state_function<T, Dim, Container>());
            state_functions.back().get()->set_base_material(make_solid_material(base_materials.back().get()));
            yield_func->state_function(0, state_functions.back().get());
            return yield_func;
        }else{
            throw std::runtime_error("general_material::make_yield_function(): no matching yield function");
        }
    }

    template<typename Json>
    constexpr inline std::unique_ptr<propagation_law_base<T>> make_propagation_law(Json const& json){
        const auto& parameter{json["parameter"]};
        if(json["name"] == "linear"){
            return make_propagation_function_linear<T>(parameter["K"]);
        }else if(json["name"] == "power_law"){
            return make_propagation_function_power_law<T>(parameter["K"], parameter["m"]);
        }else if(json["name"] == "exponential_law"){
            return make_propagation_function_exponential_law<T>(parameter["K"], parameter["m"]);
        }else if(json["name"] == "damage_exponential1"){
            return make_propagation_function_exponential1_law<T>(parameter["eps0"], parameter["b"]);
        }else if(json["name"] == "damage_exponential2"){
            if(parameter.contains("constant")){
                return make_propagation_function_exponential2_law<T>(parameter["eps0"], parameter["epsf"], parameter["constant"]);
            }else{
                return make_propagation_function_exponential2_law<T>(parameter["eps0"], parameter["epsf"], 1.0);
            }
        }else if(json["name"] == "damage_exponential3"){
            std::cout<< "damage_exponential3"<<std::endl;
            return make_propagation_function_exponential3_law<T>(parameter["kappa0"], parameter["beta"], parameter["alpha"]);
        }else{
            throw std::runtime_error("general_material::make_propagation_law(): no matching propagation law");
        }
    }

    template<typename Json>
    constexpr inline std::unique_ptr<state_function_base<T, Dim, Container>> make_state_function(Json const& json){
        if(json["name"] == "von_mises_strain"){
            return make_von_mises_strain_state_function<T, Dim, Container>();
        }else if(json["name"] == "vector_strain_state_function"){
            if(json.contains("parameter")){
                const auto& parameter{json["parameter"]};
                tmech::tensor<T, Dim, 1> m, n;
                for(std::size_t i{0}; i<Dim; ++i){
                    m(i) = parameter["m"][i];
                    n(i) = parameter["n"][i];
                }
                return make_vector_strain_state_function<T, Dim, Container>(m,n);
            }else{
                return make_vector_strain_state_function<T, Dim, Container>();
            }
        }else if(json["name"] == "strain_energy_state_function"){
            const T E{json["parameter"]["E"]};
            return make_strain_energy_strain_state_function<T, Dim, Container>(E);
        }else if(json["name"] == "von_mises"){
            return make_von_mises_state_function<T, Dim, Container>();
        }else{
            throw std::runtime_error("general_material::make_state_function(): no matching state function");
        }
    }

    template<typename Json>
    constexpr inline std::unique_ptr<eshelby_tensor_base<T, Dim>> make_eshelby_tensor(Json const& json){
        const std::string type{json["type"]};
        if(type == "long_fibre_x1"){
            auto S{make_eshelby_long_fibre_x1<T, Dim>()};
            S.get()->set_parameter(json["parameter"]["nu"]);
            return S;
        }

        if(type == "long_fibre_x2"){
            auto S{make_eshelby_long_fibre_x2<T, Dim>()};
            S.get()->set_parameter(json["parameter"]["nu"]);
            return S;
        }

        if(type == "long_fibre_x3"){
            auto S{make_eshelby_long_fibre_x3<T, Dim>()};
            S.get()->set_parameter(json["parameter"]["nu"]);
            return S;
        }

        if(type == "long_fibre"){
            auto S{make_eshelby_long_fibre<T, Dim>()};
            S.get()->set_parameter(json["parameter"]["nu"]);
            if(json["parameter"].contains("direction_vector")){
                std::vector<T> dir = json["parameter"]["direction_vector"].template get<std::vector<T>>();
                if(dir.size() != 3){
                    throw std::runtime_error("general_material::make_eshelby_tensor(): drection components != 3 long_fibre");
                }
                S.get()->direction() = tmech::adaptor<T, 3, 1, tmech::full<Dim>>(&dir[0]);
            }
            return S;
        }

        if(type == "sphere"){
            auto S{make_eshelby_solid_sphere<T, Dim>()};
            S.get()->set_parameter(json["parameter"]["nu"]);
            return S;
        }

        if(type == "ellipsoid"){
            auto S{make_eshelby_tensor_solid_cylindrical_ellipsoid<T, Dim>()};
            S.get()->set_parameter(json["parameter"]["nu"], json["parameter"]["aspec_ratio"]);
            if(json["parameter"].contains("direction_vector")){
                std::vector<T> dir = json["parameter"]["direction_vector"].template get<std::vector<T>>();
                S.get()->direction() = tmech::adaptor<T, 3, 1, tmech::full<Dim>>(&dir[0]);
            }
            return S;
        }

        if(type == "general"){
            auto S{make_eshelby_tensor_solid_cylindrical_ellipsoid<T, Dim>()};
            S.get()->set_parameter(json["parameter"]["nu"], json["parameter"]["aspec_ratio"]);
            if(json["parameter"].contains("direction_vector")){
                std::vector<T> dir = json["parameter"]["direction_vector"].template get<std::vector<T>>();
                S.get()->direction() = tmech::adaptor<T, 3, 1, tmech::full<Dim>>(&dir[0]);
            }
            return S;
        }

        //finite strain
        if(type == "finite_strain_general"){
            auto S{std::make_unique<eshelby_tensor_finite_strain_base<T,Dim,Container>>()};
            S.get()->set_parameter(json["parameter"]["a"], json["parameter"]["b"], json["parameter"]["c"]);
            return S;
        }

        //conductivity
        if(type == "conductivity_sphere"){
            return std::make_unique<eshelby_tensor_conductivity_sphere<T, Dim>>();
        }

        if(type == "conductivity_cylindrical_ellipsoid"){
            auto S{std::make_unique<eshelby_tensor_conductivity_cylindrical_ellipsoid<T, Dim>>()};
            S->set_parameter(json["parameter"]["a"], json["parameter"]["b"]);
            std::vector<T> dir = json["parameter"]["direction_vector"].template get<std::vector<T>>();
            S->direction() = tmech::adaptor<T, 3, 1, tmech::full<Dim>>(&dir[0]);
            return S;
        }

        if(type == "conductivity_cylinder"){
            auto S{std::make_unique<eshelby_tensor_conductivity_cylinder<T, Dim>>()};
            std::vector<T> dir = json["parameter"]["direction_vector"].template get<std::vector<T>>();
            S->direction() = tmech::adaptor<T, 3, 1, tmech::full<Dim>>(&dir[0]);
            return S;
        }

        if(type == "conductivity_general"){
            auto S{std::make_unique<eshelby_tensor_conductivity_general<T, Dim>>()};
            S->set_parameter(json["parameter"]["a"]);
            S->set_parameter(json["parameter"]["b"]);
            S->set_parameter(json["parameter"]["c"]);
            return S;
        }

        throw std::runtime_error("general_material::make_eshelby_tensor(): no matching eshelby tensor");
        return nullptr;
    }

    template<typename Json>
    constexpr inline std::unique_ptr<inclusion_base<T, Dim, Container>> make_inclusion(Json const& json){
        if(json["name"] == "inclusion"){
            return uvwmat::make_inclusion<T, Dim, Container>();
        }else{
            throw std::runtime_error("general_material::make_inclusion(): no matching inclusion");
        }
    }

//    template<typename Json>
//    constexpr inline void make_composite_material(composite_material_base<T, Dim, Container> * mat, Json const& json){
//        //inclusions
//        for(const auto& inclusion : json["inclusions"]){
//            inclusions.emplace_back(make_inclusion(inclusion));
//            base_materials.emplace_back(make_material_detail(inclusion["material"]));
//            eshelby_tensors.emplace_back(make_eshelby_tensor(inclusion["eshelby_tensor"]));

//            inclusions.back().get()->set_volume_fraction(inclusion["volume_fraction"]);
//            inclusions.back().get()->set_material(base_materials.back().get());
//            inclusions.back().get()->set_eshelby_tensor(eshelby_tensors.back().get());

//            mat->push_back(inclusions.back().get());
//        }

//        //matrix material
//        base_materials.emplace_back(make_material_detail(json["matrix"]["material"]));
//        mat->set_matrix_material(base_materials.back().get());
//        mat->set_matrix_volume_fraction(json["matrix"]["volume_fraction"]);
//    }

    template<typename Json>
    constexpr inline void make_mean_field_composite_material(mean_field_composite_base<T, Dim, Container> * mat, Json const& json){
        //inclusions
        std::vector<std::size_t> eshelby_idx;
        for(const auto& inclusion : json["inclusions"]){
            inclusions.emplace_back(make_inclusion(inclusion));
            //save current index, because in make_material_detail
            //new inclusions can be allocated
            const auto inc_idx{inclusions.size()-1};
            base_materials.emplace_back(make_material_detail(inclusion["material"]));
            eshelby_tensors.emplace_back(make_eshelby_tensor(inclusion["eshelby_tensor"]));
            eshelby_idx.push_back(eshelby_tensors.size() - 1);
            inclusions[inc_idx].get()->volume_fraction() = inclusion["volume_fraction"];
            inclusions[inc_idx].get()->material(base_materials.back().get());
            inclusions[inc_idx].get()->eshelby_tensor(eshelby_tensors.back().get());

            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<T,Dim,Container>*>(eshelby_tensors.back().get());
            if(S){
                S->set_base_material(base_materials.back().get());
            }

            mat->push_back(inclusions[inc_idx].get());
        }

        //matrix material
        base_materials.emplace_back(make_material_detail(json["matrix"]["material"]));
        mat->matrix_material(base_materials.back().get());
        mat->matrix_volume_fraction() = json["matrix"]["volume_fraction"];

//        //set matrix material in eshelby finite strain nonlinear
//        for(const auto i : eshelby_idx){
//            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<T,Dim,Container>*>(eshelby_tensors[i].get());
//            if(S){
//                S->set_base_material(base_materials.back().get());
//            }
//        }
    }

    constexpr inline auto arbitrary_orthogonal(tmech::tensor<value_type,3,1>const& vec)
    {
        const bool b0 = ((vec(0) <  vec(1)) && (vec(0) <  vec(2)));
        const bool b1 = ((vec(1) <= vec(0)) && (vec(1) <  vec(2)));
        const bool b2 = ((vec(2) <= vec(0)) && (vec(2) <= vec(1)));

        return tmech::cross(vec, tmech::tensor<value_type,3,1>(value_type(b0), value_type(b1), value_type(b2)));
    }


    constexpr inline auto rotation_yaw_pitch_roll(value_type const __alpha, value_type const __beta, value_type const __gamma){
        const auto sin_a{std::sin(__alpha)};
        const auto sin_b{std::sin(__beta)};
        const auto sin_g{std::sin(__gamma)};
        const auto cos_a{std::cos(__alpha)};
        const auto cos_b{std::cos(__beta)};
        const auto cos_g{std::cos(__gamma)};

        return tmech::tensor<value_type, 3, 2>{cos_b*cos_g, sin_a*sin_b*cos_g-cos_a*sin_g, cos_a*sin_b*cos_g+sin_a*sin_g,
                    cos_b*sin_g, sin_a*sin_b*sin_g+cos_a*cos_g, cos_a*sin_b*sin_g-sin_a*cos_g,
                    -sin_b,      sin_a*cos_b,                   cos_a*cos_b};
    }

    constexpr inline auto rotation_euler(value_type const __alpha, value_type const __beta, value_type const __gamma){
        return rotation_yaw_pitch_roll(__gamma, __beta, __alpha);
    }

    constexpr inline auto rotation_direction(tmech::tensor<value_type, 3, 1> const& __direction){
        //https://gamedev.stackexchange.com/questions/45298/convert-orientation-vec3-to-a-rotation-matrix
        /* Find cosφ and sinφ */
        const auto c1 = std::sqrt(__direction(0)*__direction(0) + __direction(1)*__direction(1));
        const auto s1 = __direction(2);
        /* Find cosθ and sinθ; if gimbal lock, choose (1,0) arbitrarily */
        const auto c2 = (c1 ? __direction(0) / c1 : 1.0);
        const auto s2 = (c1 ? __direction(1) / c1 : 0.0);

        return tmech::tensor<value_type, 3, 2>{__direction(0), -s2, -s1*c2,
                    __direction(1),  c2, -s1*s2,
                    __direction(2),   0,     c1};
    }

    std::unique_ptr<material_base<T, Dim, Container>>           material;
    vector<material_base<T, Dim, Container>>                    base_materials;
    vector<propagation_law_base<T>>                             propagation_laws;
    vector<yield_function_base<T, Dim ,Container>>              yield_functions;
    vector<state_function_base<T, Dim, Container>>              state_functions;
    vector<inclusion_base<T, Dim ,Container>>                   inclusions;
    vector<eshelby_tensor_base<T, Dim>>                         eshelby_tensors;
    vector<rule_of_mixture_kernal_base<T, Dim, Container>>      r_o_m_kernals;
    vector<mean_field_composite_kernal_base<T, Dim, Container>> mean_field_composite_kernals;
    vector<kinetic_model_thermoset_base<T>>                     curing_functions;
};

#endif // GENERAL_MATERIAL_BONES_H
