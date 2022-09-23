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
    json["name"] = "small_strain_mean_field_composite";
    json["kernal"] = "scs";
    json["inclusions"]["inclusion1"]["material"]["name"] = "linear_elasticity";
    json["inclusions"]["inclusion1"]["material"]["parameter"]["E"] = 2100;
    json["inclusions"]["inclusion1"]["material"]["parameter"]["nu"] = 0.3;
    json["inclusions"]["inclusion1"]["name"] = "inclusion";
    json["inclusions"]["inclusion1"]["volume_fraction"] = 0.5;
    json["inclusions"]["inclusion1"]["eshelby_tensor"]["type"] = "sphere";
    json["inclusions"]["inclusion1"]["eshelby_tensor"]["parameter"]["nu"] = 0.45;
    json["matrix"]["material"]["name"] = "linear_elasticity";
    json["matrix"]["material"]["parameter"]["E"] = 1800;
    json["matrix"]["material"]["parameter"]["nu"] = 0.45;
    json["matrix"]["volume_fraction"] = 0.5;

    uvwmat::general_material<double, Dim, std::vector<double>> general_material;
    general_material.make_material(json);
    auto* composite = uvwmat::make_composite_material(general_material.get());
    auto* material  = uvwmat::make_small_strain_material(general_material.get());

    const double E_i{json["inclusions"]["inclusion1"]["material"]["parameter"]["E"]};
    const double nu_i{json["inclusions"]["inclusion1"]["material"]["parameter"]["nu"]};
    std::cout<<"Compression modulus inclusion: "<<uvwmat::K_E_nu(E_i, nu_i)<<std::endl;
    std::cout<<"Shear modulus inclusion:       "<<uvwmat::mu_E_nu(E_i, nu_i)<<std::endl;

    const double E_m{json["matrix"]["material"]["parameter"]["E"]};
    const double nu_m{json["matrix"]["material"]["parameter"]["nu"]};
    std::cout<<"Compression modulus matrix:    "<<uvwmat::K_E_nu(E_m, nu_m)<<std::endl;
    std::cout<<"Shear modulus matrix:          "<<uvwmat::mu_E_nu(E_m, nu_m)<<std::endl;

    //plot data
    double max_steps{100};
    std::vector<value_type> c, K_eff, G_eff;
    c.reserve(max_steps);
    K_eff.reserve(max_steps);
    G_eff.reserve(max_steps);
    for(double c_m{1./max_steps}; c_m<=1.; c_m += 1./max_steps){
        c.push_back(c_m);
        composite->volume_fraction(1) = c_m;
        composite->volume_fraction(0) = 1.0 - c_m;
        material->init();
        const auto [K_MT, mu_MT]{uvwmat::extract_elasticity_parameter_K_mu(material->get_tangent_tensor())};
        K_eff.push_back(K_MT);
        G_eff.push_back(material->get_tangent_tensor()(0,1,0,1));
    }

    matplot::plot(c, K_eff);
    matplot::show();

    matplot::plot(c, G_eff);
    matplot::show();
    return 0;
}
