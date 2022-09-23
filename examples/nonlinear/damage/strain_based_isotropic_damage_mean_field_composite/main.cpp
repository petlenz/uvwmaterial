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

    json["name"] = "small_strain_mean_field_composite_damage_nonlocal";
    json["inclusions"]["inclusion1"]["material"]["name"] = "small_strain_isotropic_nonlocal_damage";
    json["inclusions"]["inclusion1"]["material"]["base_material"]["name"] = "linear_elasticity";
    json["inclusions"]["inclusion1"]["material"]["base_material"]["parameter"]["E"] = 1800;
    json["inclusions"]["inclusion1"]["material"]["base_material"]["parameter"]["nu"] = 0.45;
    json["inclusions"]["inclusion1"]["material"]["state_function"]["name"] = "von_mises_strain";
    json["inclusions"]["inclusion1"]["material"]["yield_function"]["name"] = "strain_based_damage";
    json["inclusions"]["inclusion1"]["material"]["yield_function"]["parameter"]["crit_val"] = 0.01;
    json["inclusions"]["inclusion1"]["material"]["propagation_law"]["name"] = "damage_exponential2";
    json["inclusions"]["inclusion1"]["material"]["propagation_law"]["parameter"]["eps0"] = 0.01;
    json["inclusions"]["inclusion1"]["material"]["propagation_law"]["parameter"]["epsf"] = 0.02;
    json["inclusions"]["inclusion1"]["material"]["average_type"] = "damage_variable";
    json["inclusions"]["inclusion1"]["volume_fraction"] = 0.5;
    json["inclusions"]["inclusion1"]["name"] = "inclusion";
    json["inclusions"]["inclusion1"]["eshelby_tensor"]["type"] = "sphere";
    json["inclusions"]["inclusion1"]["eshelby_tensor"]["parameter"]["nu"] = 0.45;

    json["matrix"]["material"]["name"] = "linear_elasticity";
    json["matrix"]["material"]["parameter"]["E"] = 1800;
    json["matrix"]["material"]["parameter"]["nu"] = 0.45;
    json["matrix"]["volume_fraction"] = 0.5;

    std::cout<<json<<std::endl;


    for(const auto kernal : {"dilute", "mori_tanaka", "scs"}){
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
            small_strain->strain_tensor() = eps*i;

            //update tangent and stresses
            small_strain->update();

            //get constant reference of stress tensor
            const auto& sigma{small_strain->stress_tensor()};
            const auto& sigma_i{composite->material(0)->stress_tensor()};
            const auto& sigma_m{composite->material(1)->stress_tensor()};

            //determine von mises stress
            von_Mises.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma),tmech::dev(sigma))));
            von_Mises_m.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma_m),tmech::dev(sigma_m))));
            von_Mises_i.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma_i),tmech::dev(sigma_i))));

            //history
            D_i.push_back(history->history()[0]);
        }

//        matplot::plot(load_steps, von_Mises, "-o -black",
//                      load_steps, von_Mises_i, "-x -blue",
//                      load_steps, von_Mises_m, "-* -red");
//        matplot::title(kernal);
//        matplot::show();

//        matplot::plot(load_steps, D_i, "-o -black");
//        matplot::title(kernal);
//        matplot::show();
    }

    return 0;
}
