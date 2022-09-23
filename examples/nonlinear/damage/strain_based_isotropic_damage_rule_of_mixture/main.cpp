#include <iostream>
#include <matplot/matplot.h>
#include <tmech/tmech.h>
#include "../../../../include/uvwmat.h"
#include "../../../../external_libraries/nlohmann/json.hpp"


int main()
{
    constexpr std::size_t Dim{3};
    using value_type = double;
    nlohmann::json json;

    json["name"] = "small_strain_rule_of_mixture_history";
    //json["kernal"] = "reuss";
    json["materials"]["material1"]["name"] = "small_strain_isotropic_damage";
    json["materials"]["material1"]["base_material"]["name"] = "linear_elasticity";
    json["materials"]["material1"]["base_material"]["parameter"]["E"] = 1800;
    json["materials"]["material1"]["base_material"]["parameter"]["nu"] = 0.45;
    json["materials"]["material1"]["state_function"]["name"] = "von_mises_strain";
    json["materials"]["material1"]["yield_function"]["name"] = "strain_based_damage";
    json["materials"]["material1"]["yield_function"]["parameter"]["crit_val"] = 0.01;
    json["materials"]["material1"]["propagation_law"]["name"] = "damage_exponential2";
    json["materials"]["material1"]["propagation_law"]["parameter"]["eps0"] = 0.01;
    json["materials"]["material1"]["propagation_law"]["parameter"]["epsf"] = 0.02;
    json["materials"]["material1"]["volume_fraction"] = 0.5;

    json["materials"]["material2"]["name"] = "linear_elasticity";
    json["materials"]["material2"]["parameter"]["E"] = 1800;
    json["materials"]["material2"]["parameter"]["nu"] = 0.45;
    json["materials"]["material2"]["volume_fraction"] = 0.5;

    for(const auto kernal : {"reuss", "voigt", "vrh"}){
        json["kernal"] = kernal;
        uvwmat::general_material<double, Dim, std::vector<double>> general_material;
        general_material.make_material(json);
        general_material.get()->init();
        auto small_strain = uvwmat::make_small_strain_material(general_material.get());
        auto composite = uvwmat::make_composite_material(general_material.get());
        auto history = uvwmat::make_history(general_material.get());

        std::vector<double> load_steps, von_Mises, von_Mises_m, von_Mises_i, D_i;
        tmech::tensor<double, Dim, 2> eps;
        eps(0,0) = 0.05;

        for(double i{0}; i<1.0; i += 0.002){
            //for plotting
            load_steps.push_back(i);

            //set strain tensor
            small_strain->set_strain_tensor(eps*i);

            //update tangent and stresses
            small_strain->update();

            //get constant reference of stress tensor
            const auto& sigma{small_strain->get_stress_tensor()};
            const auto& sigma_i{composite->material(0)->get_stress_tensor()};
            const auto& sigma_m{composite->material(1)->get_stress_tensor()};

            //determine von mises stress
            von_Mises.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma),tmech::dev(sigma))));
            von_Mises_m.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma_m),tmech::dev(sigma_m))));
            von_Mises_i.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma_i),tmech::dev(sigma_i))));

            //history
            D_i.push_back(history->get_history()[0]);
        }

        matplot::plot(load_steps, von_Mises, "-o -black",
                      load_steps, von_Mises_i, "-x -blue",
                      load_steps, von_Mises_m, "-* -red");
        matplot::title(kernal);
        matplot::show();

        matplot::plot(load_steps, D_i, "-o -black");
        matplot::title(kernal);
        matplot::show();
    }

    return 0;
}
