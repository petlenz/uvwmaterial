#include <iostream>
#include <matplot/matplot.h>
#include <tmech/tmech.h>
#include <boost/container/flat_map.hpp>
#include <boost/math/special_functions.hpp>
#include <uvwmaterial/uvwmaterial.h>
#include <nlohmann/json.hpp>

using namespace std;

static constexpr std::size_t Dim{3};
using value_type = double;
using tensor2    = tmech::tensor<value_type, Dim, 2>;
using tensor4    = tmech::tensor<value_type, Dim, 4>;

int main()
{
    tensor2 eps;
    eps(0,0) = 1;
    
    nlohmann::json json;
    json["name"] = "linear_elasticity";
    json["parameter"]["E"] = 2100;
    json["parameter"]["nu"] = 0.3;

    uvwmat::general_material<value_type, Dim, std::vector<value_type>> mat;
    mat.make_material(json);

    mat.get()->init();
    auto small_strain = uvwmat::make_small_strain_material(mat.get());
    std::vector<double> load_steps, von_Mises;

    for(double i{0}; i<1.0; i += 0.1){
        //for plotting
        load_steps.push_back(i);

        //set strain tensor
        small_strain->strain_tensor() = eps*i;

        //update tangent and stresses
        small_strain->update();

        //get constant reference of stress tensor
        const auto& sigma{small_strain->stress_tensor()};

        //determine von mises stress
        von_Mises.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma),tmech::dev(sigma))));
    }

    matplot::plot(load_steps, von_Mises, "-o -black");
    matplot::show();


    return 0;
}
