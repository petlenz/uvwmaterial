#include <iostream>
#include <matplot/matplot.h>
#include <tmech/tmech.h>
#include "../../../../include/uvwmat.h"
#include "../../../../external_libraries/nlohmann/json.hpp"

//composite
//--> virtual number_of_materials()
//--> virtual material(number)
//--> virtual volume_fraction(number)

//history_composite_base : history_base
//--> pointer to composite mateiral

//nonlocal_composite_base : nonlocal_base
//--> pointer to composite mateiral


//rule_of_mixture : composite
//--> contains data
//--> usable for small and large strains

//small_strain_rule_of_mixture : rule_of_mixture, small_strain_material
//--> update()
//--> update_stress()
//--> update_tangent()



//mean_field_composite : composite
//--> contains data
//--> usable for small and large strains




int main()
{
    constexpr std::size_t Dim{3};
    using value_type = double;
    nlohmann::json json;
    json["name"] = "small_strain_isotropic_damage";
    json["base_material"]["name"] = "linear_elasticity";
    json["base_material"]["parameter"]["E"] = 1800;
    json["base_material"]["parameter"]["nu"] = 0.45;
    json["state_function"]["name"] = "von_mises_strain";
    json["yield_function"]["name"] = "strain_based_damage";
    json["yield_function"]["parameter"]["crit_val"] = 0.01;
    json["propagation_law"]["name"] = "damage_exponential2";
    json["propagation_law"]["parameter"]["eps0"] = 0.01;
    json["propagation_law"]["parameter"]["epsf"] = 0.02;

    uvwmat::general_material<double, Dim, std::vector<double>> general_material;
    general_material.make_material(json);
    general_material.get()->init();
    auto small_strain = uvwmat::make_small_strain_material(general_material.get());
    std::vector<double> load_steps, von_Mises;
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

        //determine von mises stress
        von_Mises.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma),tmech::dev(sigma))));
    }

    matplot::plot(load_steps, von_Mises, "-o -black");
    matplot::show();

    return 0;
}
