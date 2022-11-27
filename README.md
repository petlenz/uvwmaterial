## Introduction

`uvwmaterial` is a C++ library for the numerical study of the physics of continuous materials.

`uvwmaterial` provides

- tba


`uvwmaterial` requires a modern C++ compiler supporting C++17. The following C++
compilers are supported:

- On Unix platforms, gcc 5 or a recent version of Clang
- On Windows platforms not tested

## Install from sources

tba

## Documentation

Current supported materials:

- Linear elasticity

For more information on using `uvwmaterial`, check out the documentation

https://uvwmaterial.readthedocs.io

## Issues and enhancements

We use the [GitHub issue
tracker](https://github.com/petlenz/uvwmaterial/issues) for all bug/issue
reports and proposals for enhancements.

## Usage

tba

### Basic usage

..code 
#include <iostream>
#include <matplot/matplot.h>
#include <tmech/tmech.h>
#include <boost/container/flat_map.hpp>
#include <uvwmaterial/uvwmaterial.h>
```cpp
int main()
{
    double E{210000}, nu{0.33};
    uvwmat::linear_elasticity<double, 3, std::vector<double>> linear_elasticity(E, nu);

    //or
    //linear_elasticity.set_parameter(std::vector<double>{E,nu});

    //or
    //linear_elasticity.set_parameter(E, nu);

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
        linear_elasticity.strain_tensor() = eps*i;

        //update tangent and stresses
        linear_elasticity.update();

        //get constant reference of stress tensor
        const auto& sigma{linear_elasticity.stress_tensor()};

        //determine von mises stress
        von_Mises.push_back(std::sqrt(1.5*tmech::dcontract(tmech::dev(sigma),tmech::dev(sigma))));
    }

    matplot::plot(load_steps, von_Mises, "-o -black");
    matplot::xlabel("load steps");
    matplot::ylabel("von Mises stress");
    matplot::title("Linear elasticity");
    matplot::save("linear_elasticity.png");
    //matplot::show();


    return 0;
}
```
