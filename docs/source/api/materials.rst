.. Copyright (c) 2022, Peter Lenz

   Distributed under the terms of the BSD 3-Clause License.

   The full license is in the file LICENSE, distributed with this software.


=========
Materials
=========


Nonlocal
========

Damage
------

Plasticity
----------

Linear
======

Json Linear Elasticity
----------------------

.. code::

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


Linear Elasticity
-----------------

.. code::

    #include <iostream>
    #include <matplot/matplot.h>
    #include <tmech/tmech.h>
    #include <boost/container/flat_map.hpp>
    #include <uvwmaterial/uvwmaterial.h>
    
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


Composite
=========

Nonlocal
--------

Linear
------

Nonlinear
---------

Nonlinear
=========

Damage
------

Plasticity
----------

