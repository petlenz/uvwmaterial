#include <iostream>
#include <matplot/matplot.h>
#include <tmech/tmech.h>
#include "../../../../../include/uvwmat.h"
#include "../../../../../external_libraries/nlohmann/json.hpp"


int main()
{
    constexpr std::size_t Dim{3};
    using value_type = double;
    nlohmann::json json;
    json["name"] = "small_strain_short_fibre_composite";
    json["integration_precision"] = 7;
    json["fibre_orientation_tensor"] = std::vector<double>{1./3., 0, 0, 0, 1./3., 0, 0, 0, 1./3.};

    json["composite"]["name"] = "small_strain_mean_field_composite";
    json["composite"]["kernal"] = "mori_tanaka";
    json["composite"]["inclusions"]["inclusion1"]["material"]["name"] = "linear_elasticity";
    json["composite"]["inclusions"]["inclusion1"]["material"]["parameter"]["E"] = 2100;
    json["composite"]["inclusions"]["inclusion1"]["material"]["parameter"]["nu"] = 0.3;
    json["composite"]["inclusions"]["inclusion1"]["name"] = "inclusion";
    json["composite"]["inclusions"]["inclusion1"]["volume_fraction"] = 0.5;
    json["composite"]["inclusions"]["inclusion1"]["eshelby_tensor"]["type"] = "ellipsoid";
    json["composite"]["inclusions"]["inclusion1"]["eshelby_tensor"]["parameter"]["nu"] = 0.45;
    json["composite"]["inclusions"]["inclusion1"]["eshelby_tensor"]["parameter"]["aspec_ratio"] = 27.;
    json["composite"]["matrix"]["material"]["name"] = "linear_elasticity";
    json["composite"]["matrix"]["material"]["parameter"]["E"] = 1800;
    json["composite"]["matrix"]["material"]["parameter"]["nu"] = 0.45;
    json["composite"]["matrix"]["volume_fraction"] = 0.5;

    uvwmat::general_material<double, Dim, std::vector<double>> general_material;
    general_material.make_material(json);
    general_material.get()->init();
    auto small_strain = uvwmat::make_small_strain_material(general_material.get());
    std::vector<double> load_steps, von_Mises;
    tmech::tensor<double, Dim, 2> eps;
    eps(0,0) = 1;

    for(double i{0}; i<1.0; i += 0.1){
        //for plotting
        load_steps.push_back(i);

        //set strain tensor
        small_strain->set_strain_tensor(eps*i);

        //update tangent and stresses
        small_strain->update();

        //get constant reference of stress tensor
        const auto& sigma{small_strain->get_stress_tensor()};

        //determine von mises stress
        von_Mises.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma),tmech::dev(sigma))));
    }

    matplot::plot(load_steps, von_Mises, "-o -black");
    matplot::show();
    return 0;
}
