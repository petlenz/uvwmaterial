#include <iostream>
#include <matplot/matplot.h>
#include <tmech/tmech.h>
#include "../../../include/uvwmat.h"

int main()
{
    double E{210000}, nu{0.33};
    uvwmat::linear_elasticity<double, 3, std::vector<double>> linear_elasticity(E, nu);

    //or
    //linear_elasticity.set_parameter(std::vector<double>{E,nu});

    //or
    //linear_elasticity.set_parameter(E, nu);

    //or
    //std::vector<double> parameter{E, nu};
    //linear_elasticity.set_parameter(parameter.begin(), parameter.end());

    //or
    //linear_elasticity.set_parameter(std::make_pair(0ul, E));
    //linear_elasticity.set_parameter(std::make_pair(1ul, nu));

    tmech::tensor<double, 3, 2> eps;
    eps(0,0) = 1;

    //init material;
    //setup elasticity tensor
    linear_elasticity.init();

    std::vector<double> load_steps, von_Mises;

    for(double i{0}; i<1.0; i += 0.1){
        //for plotting
        load_steps.push_back(i);

        //set strain tensor
        linear_elasticity.set_strain_tensor(eps*i);

        //update tangent and stresses
        linear_elasticity.update();

        //get constant reference of stress tensor
        const auto& sigma{linear_elasticity.get_stress_tensor()};

        //determine von mises stress
        von_Mises.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma),tmech::dev(sigma))));
    }

    matplot::plot(load_steps, von_Mises, "-o -black");
    matplot::show();


    return 0;
}
