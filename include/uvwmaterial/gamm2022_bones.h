#ifndef GAMM2022_BONES_H
#define GAMM2022_BONES_H

#include "material_base_bones.h"
#include "thermo_mechanical_material_base_bones.h"
#include "solid_material_base_bones.h"
#include "conductivity_material_base_bones.h"
#include "polymer_curing_functions_bones.h"

template <typename _T, std::size_t _Dim, typename _Container>
class gamm2022 :
        public thermo_mechanical_material_base<_T, _Dim, _Container>,
        public incremental_solid_material<_T, _Dim>,
        public history_material_base<_T>,
        public time_dependent_material_base<_T>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>
{
public:

    using size_type = std::size_t;

    gamm2022():
        thermo_mechanical_material_base<_T, _Dim, _Container>(),
        incremental_solid_material<_T, _Dim>(),
        history_material_base<_T>(),
        time_dependent_material_base<_T>(),
        finite_strain_solid_material_base<_T, _Dim, _Container>(this),
        _solid_mat(nullptr),
        _thermal_mat(nullptr),
        _curing_function(nullptr),
        _this_history(0)
    {}

    virtual ~gamm2022() {}

    inline virtual void init(){
        if(!this->_is_init){
            if(!_solid_mat || !_thermal_mat || !_curing_function){
                throw std::runtime_error("thermo_mechanical_material_base::init() missing data");
            }

            _solid_mat->init();
            _thermal_mat->init();

            this->_type = FINITE_STRAIN_FORMULATION::MIXED;
            //degree of cure
            _this_history += 18;

            //all history
            size_type his_size{0};

            auto thermal_his = dynamic_cast<history_material_base<_T>*>(_thermal_mat);
            if(thermal_his){
                his_size += thermal_his->size();
            }

            auto solid_his = dynamic_cast<history_material_base<_T>*>(_solid_mat);
            if(solid_his){
                his_size += solid_his->size();
            }

            history_material_base<_T>::_history.resize(_this_history+his_size);


            //            _curing_mat.resize(1);
            //            _temperature_mat.resize(4);
            //            _composite_solid_mat.resize(3);
            //            _composite_conductivity_mat.resize(3);



            get_history_local();

            this->_is_init = true;
        }
    }

    inline virtual void reinit(){
        _solid_mat->reinit();
        _thermal_mat->reinit();
    }

    inline virtual void update(){

        set_history_local();

        //update volume fractions
        const auto n{this->_parameter[0]};
        const auto cs{this->_history[0]};
        const auto cc{(1-n)*(1-cs)};
        const auto cr{n*(1-cs)};

        //this->_temperature += 273.15;


        //update volume fractions in thermal material
        // composite {inclusions, matrix} matrix -> polymer composite
        auto thermal_composite{dynamic_cast<mean_field_composite_base<_T, _Dim, _Container>*>(_thermal_mat)};
        //get matrix
        auto thermal_matrix_composite = dynamic_cast<mean_field_composite_base<_T, _Dim, _Container>*>(thermal_composite->material(thermal_composite->number_of_materials()-1));
        //get inclusion sc
        auto thermal_sc_composite = dynamic_cast<mean_field_composite_base<_T, _Dim, _Container>*>(thermal_matrix_composite->material(0));

        //inclusion {solid+curing agent} matrix resin
        thermal_matrix_composite->volume_fraction(0) = 1.0 - cr; //inclusion
        thermal_matrix_composite->volume_fraction(1) = cr; //matrix
        //inclusion solid matrix curing agent
        thermal_sc_composite->volume_fraction(0) = cs/(cc + cs); //inclusion
        thermal_sc_composite->volume_fraction(1) = cc/(cc + cs); //matrix


        //general_finite_strain_conductivity
        //conductivity_mean_field explicit
        //conductivity_finite_strain


        //update geometry thermal conductivity
        //R_3 = Rn*(cs+cc+cr)^(1/3); --> radius sphere
        //R_2 = Rn*(cs+cc)^{1/3};    --> radius resin
        //R_1 = Rn*(cs)              --> radius solid
        {
            auto& para_1 = thermal_matrix_composite->inclusion(0)->eshelby_tensor()->parameter();
            const auto scale{std::min((_T)1.0,this->_parameter[4])};
            const auto Rn{this->_parameter[4]/scale};
            para_1[0] = (cs+cc == 0 ? 1e-12 : Rn*std::pow(cs+cc, 1./3.)*scale);
            para_1[1] = para_1[0];
            para_1[2] = para_1[0];
            auto& para_2 = thermal_sc_composite->inclusion(0)->eshelby_tensor()->parameter();
            para_2[0] = (cs == 0 ? 1e-12 : Rn*std::pow(cs, 1./3.)*scale);
            para_2[1] = para_2[0];
            para_2[2] = para_2[0];
            this->_history[4] = this->_parameter[4];
            this->_history[5] = para_1[0];
            this->_history[6] = para_2[0];
            //std::cout<<cs<<" "<<cc<<" "<<cr<<std::endl;
            //std::cout<<para_1[0]<<" "<<para_2[0]<<std::endl;
            //            if(para_1[1] < para_2[1]){
            //                throw std::runtime_error("gamm2022::update() radien are equal");
            //            }
        }


        this->_history[3] = dynamic_cast<conductivity_material_base<_T, _Dim, _Container>*>(thermal_matrix_composite)->conductivity_tensor()(0,0);


        _thermal_mat->reinit(); //nötig?
        //update heat flux, conductivity, ...
        _thermal_mat->scalar_value() = this->_temperature;
        _thermal_mat->gradient_tensor() = this->thermal_gradient();
        _thermal_mat->update();

        this->_q = tmech::det(this->_Fn)*tmech::inv(this->_Fn)*_thermal_mat->flux_tensor();
        this->_k = tmech::det(this->_Fn)*tmech::inv(this->_Fn)*_thermal_mat->conductivity_tensor()*tmech::trans(this->_Fn);
        this->_q = _thermal_mat->flux_tensor();
        this->_k = _thermal_mat->conductivity_tensor();

        //this->_q.print(std::cout);
        //this->_k.print(std::cout);

        //update volume fraction in solid material
        auto solid_composite{dynamic_cast<mean_field_composite_base<_T, _Dim, _Container>*>(_solid_mat)};
        //get matrix
        auto solid_matrix_composite = dynamic_cast<mean_field_composite_base<_T, _Dim, _Container>*>(solid_composite->material(solid_composite->number_of_materials()-1));
        //get inclusion sc
        auto solid_sc_composite = dynamic_cast<mean_field_composite_base<_T, _Dim, _Container>*>(solid_matrix_composite->material(0));



        //update volume fractions
        solid_matrix_composite->volume_fraction(0) = cc + cs; //inclusion
        solid_matrix_composite->volume_fraction(1) = cr; //matrix
        solid_sc_composite->volume_fraction(0) = cs/(cc + cs); //inclusion
        solid_sc_composite->volume_fraction(1) = cc/(cc + cs); //matrix

        //        std::cout<<"Old volume fraction\n";
        //        std::cout<<"cs "<<cs<<" cc "<<cc<<" cr "<<cr<<" "<<cs+cr+cc<<std::endl;
        //        std::cout<<1.0 - cr<<" "<<cr<<std::endl;
        //        std::cout<<cs/(cc + cs)<<" "<<cc/(cc + cs)<<" "<<cc/(cc + cs) + cs/(cc + cs)<<std::endl;

        //inclusion material
        auto inc_mat{static_cast<solid_material_base<_T,_Dim,_Container>*>(solid_composite->material(0))};
        //resin material
        auto resin_mat{static_cast<solid_material_base<_T,_Dim,_Container>*>(solid_matrix_composite->material(1))};
        //curing agent material
        auto curing_agent_mat{static_cast<solid_material_base<_T,_Dim,_Container>*>(solid_sc_composite->material(1))};
        //solid material
        auto solid_sphere_mat{static_cast<solid_material_base<_T,_Dim,_Container>*>(solid_sc_composite->material(0))};


        //update temperature T^n
        //inclusion material
        auto inc_temp_mat{dynamic_cast<temperature_dependent_material_base<_T>*>(inc_mat)};
        if(inc_temp_mat){
            inc_temp_mat->temperature() = this->_temperature - this->_delta_temperature;
        }
        //resin material
        auto resin_temp_mat{dynamic_cast<temperature_dependent_material_base<_T>*>(resin_mat)};
        if(resin_temp_mat){
            resin_temp_mat->temperature() = this->_temperature - this->_delta_temperature;
        }
        //curing agent material
        auto curing_temp_agent_mat{dynamic_cast<temperature_dependent_material_base<_T>*>(curing_agent_mat)};
        if(curing_temp_agent_mat){
            curing_temp_agent_mat->temperature() = this->_temperature - this->_delta_temperature;
        }
        //solid material
        auto solid_sphere_temp_mat{dynamic_cast<temperature_dependent_material_base<_T>*>(solid_sphere_mat)};
        if(solid_sphere_temp_mat){
            solid_sphere_temp_mat->temperature() = this->_temperature - this->_delta_temperature;
        }

        //update curing z^n
        auto solid_sphere_curing_mat{dynamic_cast<degree_of_cure_dependent_material_base<_T>*>(solid_sphere_mat)};
        if(solid_sphere_curing_mat){
            solid_sphere_curing_mat->degree_of_cure() = this->_history[0];
            solid_sphere_curing_mat->delta_degree_of_cure() = this->_history[1];
        }


        //update solid composite
        dynamic_cast<incremental_solid_material<_T, _Dim>*>(_solid_mat)->dstrain_tensor() = this->_dstrain;
        auto finite_strain_mean_field = dynamic_cast<finite_strain_mean_field_composite_explicit<_T,_Dim,_Container>*>(_solid_mat);
        finite_strain_mean_field->set_local_history();
        finite_strain_mean_field->update_strain();
        finite_strain_mean_field->update_material();
        finite_strain_mean_field->update_macroscopic_properties();



        //this->_C.print(std::cout);
        //this->_stress.print(std::cout);

        //volume fraction inclusion
        const auto ci{solid_composite->volume_fraction(0)};
        //volume fraction matrix
        const auto cm{solid_composite->volume_fraction(1)};

        //heat capacity
        this->_cd = cm*this->_parameter[1] + ci*this->_parameter[5];

        //dichte
        this->_rho = cm*this->_parameter[2] + ci*this->_parameter[6];

        this->_latent_heat = -cm*this->_parameter[3]*this->_parameter[2]*this->_history[2];


        //update macro stress and tangent
        this->_C = _solid_mat->tangent_tensor();
        this->_stress = _solid_mat->stress_tensor();

        //this->_history[7] = tmech::dcontract(this->_stress, this->_dstrain)/this->_delta_time;
        //std::cout<<"power "<<this->_history[7]<<std::endl;
        //std::cout<<"latent heat "<<this->_latent_heat<<std::endl;


        //von mises inclusion
        {
            auto mat = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(inc_mat);
            const auto sig_dev = tmech::dev(mat->spatial_stress_tensor());
            this->_history[7] = std::sqrt(1.5*tmech::dcontract(sig_dev, sig_dev));
        }

        //von mieses matrix
        {
            auto mat = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(solid_matrix_composite);
            const auto sig_dev = tmech::dev(mat->spatial_stress_tensor());
            this->_history[8] = std::sqrt(1.5*tmech::dcontract(sig_dev, sig_dev));
        }

        //update curing
        _curing_function->set_temperature(this->_temperature+273.15);
        _curing_function->set_history(this->_history[0]);
        _curing_function->set_time_increment(this->_delta_time);
        _curing_function->update();
        const auto old_curing{this->_history[0]};
        this->_history[0] = _curing_function->get_degree_of_cure();
        this->_history[1] = _curing_function->get_degree_of_cure() - old_curing;
        this->_history[2] = (this->_history[0]-old_curing)/this->_delta_time;
        this->_history[3] = cm*this->_history[1]*this->_parameter[2]*this->_parameter[3];



        //update solid concentration tensors
        {
            const auto cs{this->_history[0]};
            const auto cc{(1-n)*(1-cs)};
            const auto cr{n*(1-cs)};

            //update volume fractions
            solid_matrix_composite->volume_fraction(0) = cc + cs; //inclusion
            solid_matrix_composite->volume_fraction(1) = cr; //matrix
            solid_sc_composite->volume_fraction(0) = cs/(cc + cs); //inclusion
            solid_sc_composite->volume_fraction(1) = cc/(cc + cs); //matrix

            //            std::cout<<"New volume fraction\n";
            //            std::cout<<"cs "<<cs<<" cc "<<cc<<" cr "<<cr<<std::endl;
            //            std::cout<<1.0 - cr<<" "<<cr<<std::endl;
            //            std::cout<<cs/(cc + cs)<<" "<<cc/(cc + cs)<<" "<<cc/(cc + cs) + cs/(cc + cs)<<std::endl;

            auto& para_1 = solid_matrix_composite->inclusion(0)->eshelby_tensor()->parameter();
            const auto scale{std::min((_T)1.0,this->_parameter[4])};
            const auto Rn{this->_parameter[4]/scale};
            para_1[0] = (cs+cc == 0 ? 1e-12 : Rn*std::pow(cs+cc, 1./3.)*scale);
            para_1[1] = para_1[0];
            para_1[2] = para_1[0];
            auto& para_2 = solid_sc_composite->inclusion(0)->eshelby_tensor()->parameter();
            para_2[0] = (cs == 0 ? 1e-8 : Rn*std::pow(cs, 1./3.)*scale);
            para_2[1] = para_2[0];
            para_2[2] = para_2[0];

            //            std::cout<<"cs "<<cs<<" cc "<<cc<<" cr "<<cr<<std::endl;
            //            std::cout<<"Rr "<<para_1[1]<<" Rs "<<para_2[1]<<std::endl;
            //            if(para_1[1] < para_2[1]){
            //                throw std::runtime_error("gamm2022::update() radien are equal");
            //            }
        }
        finite_strain_mean_field->update_strain_concentration_tensor();
        finite_strain_mean_field->get_local_history();

        {
            {
                //solid + curing agent as inclusion
                auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(solid_matrix_composite->inclusion(0)->eshelby_tensor());
                //const auto& R = S->rotation_tensor();
                //auto& para_3 = S->updated_parameter();
                tmech::adaptor<_T,_Dim,1,tmech::full<_Dim>> his(&this->_history[9]);
                //his = tmech::abs(tmech::adaptor<_T,_Dim,1,tmech::full<_Dim>>(&para_3[0]));
                his = S->geometry();
            }
            {
                //solid as inclusion
                auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(solid_sc_composite->inclusion(0)->eshelby_tensor());
                //const auto& R = S->rotation_tensor();
                //auto& para_3 = S->updated_parameter();
                tmech::adaptor<_T,_Dim,1,tmech::full<_Dim>> his(&this->_history[12]);
                his = S->geometry();
                //his = tmech::abs(tmech::adaptor<_T,_Dim,1,tmech::full<_Dim>>(&para_3[0]));
            }
            {
                //inclusion
                auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(solid_composite->inclusion(0)->eshelby_tensor());
                //const auto& R = S->rotation_tensor();
                //auto& para_3 = S->updated_parameter();
                tmech::adaptor<_T,_Dim,1,tmech::full<_Dim>> his(&this->_history[15]);
                his = S->geometry();//tmech::abs(tmech::adaptor<_T,_Dim,1,tmech::full<_Dim>>(&para_3[0]));
            }
        }
        get_history_local();
    }

    inline virtual void update_stress(){

    }

    inline virtual void update_tangent(){

    }

    inline virtual void update_heat_flux(){

    }

    inline virtual void update_conductivity(){

    }

    inline virtual thermo_mechanical_material_base<_T, _Dim, _Container>* base_material() override{
        return this;
    }

    constexpr inline auto set_solid_material(material_base<_T,_Dim,_Container>* __solid_mat){
        _solid_mat = static_cast<solid_material_base<_T,_Dim,_Container>*>(__solid_mat);
    }

    constexpr inline auto set_thermo_material(material_base<_T,_Dim,_Container>* __thermal_mat){
        _thermal_mat = static_cast<conductivity_material_base<_T,_Dim,_Container>*>(__thermal_mat);
    }

    constexpr inline auto set_curing_function(kinetic_model_thermoset_base<_T>* __curing_function){
        _curing_function = __curing_function;
    }


    constexpr inline auto get_history_local(){
        size_type iter{0};

        auto thermal_his = dynamic_cast<history_material_base<_T>*>(_thermal_mat);
        if(thermal_his){
            for(const auto& his : thermal_his->history()){
                this->_history[_this_history+iter++] = his;
            }
        }

        auto solid_his = dynamic_cast<history_material_base<_T>*>(_solid_mat);
        if(solid_his){
            for(const auto& his : solid_his->history()){
                this->_history[_this_history+iter++] = his;
            }
        }
    }

    constexpr inline auto set_history_local(){
        size_type iter{0};

        auto thermal_his = dynamic_cast<history_material_base<_T>*>(_thermal_mat);
        if(thermal_his){
            for(auto& his : thermal_his->history()){
                his = this->_history[_this_history+iter++];
            }
        }

        auto solid_his = dynamic_cast<history_material_base<_T>*>(_solid_mat);
        if(solid_his){
            for(auto& his : solid_his->history()){
                his = this->_history[_this_history+iter++];
            }
        }

    }

protected:
    solid_material_base<_T,_Dim,_Container>* _solid_mat;
    conductivity_material_base<_T,_Dim,_Container>* _thermal_mat;
    kinetic_model_thermoset_base<_T>* _curing_function;
    std::size_t _this_history;

    std::vector<degree_of_cure_dependent_material_base<_T>> _curing_mat;
    std::vector<temperature_dependent_material_base<_T>> _temperature_mat;
    std::vector<composite_material_solid_base<_T, _Dim, _Container>> _composite_solid_mat;
    std::vector<composite_material_conductivity_base<_T, _Dim, _Container>> _composite_conductivity_mat;
};

#endif // GAMM2022_BONES_H
