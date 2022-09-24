/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
﻿#ifndef FINITE_STRAIN_MEAN_FIELD_COMPOSITE_BONES_H
#define FINITE_STRAIN_MEAN_FIELD_COMPOSITE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class finite_strain_mean_field_composite :
        public mean_field_composite_solid_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>,
        public incremental_solid_material<_T, _Dim>,
        public composite_history_base<_T, _Dim, _Container>,
        public solid_material_base<_T, _Dim, _Container>
{
    //Mixed <vol_frac, solid_material, strain_ct>
    struct mixed_formulation_wrapper
    {
        _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        finite_strain_solid_material_base<_T, _Dim, _Container>* _finite_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
    //Reference <vol_frac, solid_material, finite_material, strain_ct>
    struct reference_formulation_wrapper
    {
        _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        finite_strain_solid_material_base<_T, _Dim, _Container>* _finite_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
    using adaptor = tmech::adaptor<_T,_Dim, 2, tmech::full<_Dim>>;
    using tensor2 = tmech::tensor<_T, _Dim, 2>;
public:
    using value_type = _T;
    using size_type  = std::size_t;
    using history_base = composite_history_base<_T, _Dim, _Container>;

    finite_strain_mean_field_composite():
        mean_field_composite_solid_base<_T, _Dim, _Container>(),
        finite_strain_solid_material_base<_T, _Dim, _Container>(this),
        incremental_solid_material<_T, _Dim>(),
        composite_history_base<_T, _Dim, _Container>(this),
        _finite_strain_materials(),
        _solid_materials(),
        _mixed_formulations(),
        _reference_formulations(),
        _local_history_size(0)
    {}

    virtual ~finite_strain_mean_field_composite() {}

    inline virtual void init() override {
        if(!solid_material_base<_T, _Dim, _Container>::_is_init || !mean_field_composite_solid_base<_T, _Dim, _Container>::_is_init){
            this->_matrix_material->init();
            for(auto inclusion : this->_inclusions){
                inclusion->init();
            }

            finite_strain_solid_material_base<_T, _Dim, _Container>::_type = FINITE_STRAIN_FORMULATION::MIXED;

            mean_field_composite_solid_base<_T, _Dim, _Container>::init();


            //setup finite strain material pointer
            _finite_strain_materials.resize(this->_inclusions.size()+1);
            _finite_strain_materials.back() = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _finite_strain_materials[i] = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }


            //setup solid material pointer
            _solid_materials.resize(this->_inclusions.size()+1);
            _solid_materials.back() = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _solid_materials[i] = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }

            //setup mixed formulation wrappers
            //and init F_i = I
            for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
                _finite_strain_materials[i]->deformation_tensor() = tmech::eye<value_type,_Dim,2>();
                auto& vol_frac{(i == _finite_strain_materials.size() - 1 ? this->_cm : this->_inclusions[i]->volume_fraction())};
                switch (_finite_strain_materials[i]->formulation_type()) {
                case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
                    _reference_formulations.push_back({vol_frac, _solid_materials[i], _finite_strain_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                case FINITE_STRAIN_FORMULATION::MIXED:
                    _mixed_formulations.push_back({vol_frac, _solid_materials[i], _finite_strain_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                }
            }

            //init history
            //get local history size
            history_base::init();
            _local_history_size = history_base::size();
            auto his_size{history_base::size()};
            //entries of each deformation tensor
            his_size += _Dim*_Dim*_solid_materials.size();




            //initial values
            //vol fraction at reference solution
            auto init_val{_solid_materials.size()};
            //for each finite_strain_eshelby_tensor
            //-->  geometry at tn
            for(auto inclusion : this->_inclusions){
                auto S = inclusion->eshelby_tensor();
                if(dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(S)){
                    init_val += 3;//geometry + nu
                }
            }

            //get init volume fraction, geometry
            _init_values.resize(init_val);
            get_initial_values();

            //resize history
            history_base::resize(his_size+init_val);//init_val only for illustration

            get_local_history();

            update_macro_tangent();

            solid_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    inline virtual void reinit() override {
        this->_matrix_material->reinit();
        for(auto inclusion : this->_inclusions){
            inclusion->reinit();
        }
        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();
        update_macro_tangent();
    }

    inline virtual void update_strain()override{
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_strain(): class is not initialized");
        }

        //const auto& F_macro{this->deformation_tensor() + this->_dstrain};
        for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
            _finite_strain_materials[i]->deformation_tensor() += tmech::dcontract(this->_strain_concentration_tensors[i], this->_dstrain);
            //_finite_strain_materials[i]->deformation_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], F_macro);
        }
    }

    inline virtual void update() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update(): is not initialized");
        }


        std::vector<tensor2> R(_solid_materials.size());
        std::vector<tensor2> dFi(_solid_materials.size());

        size_type iter{0}, max_iter{static_cast<size_type>(this->_parameter[1])};
        value_type tol{this->_parameter[0]};

        set_local_history();
        set_initial_values();

        for(auto solid : _solid_materials){
            solid->reinit();
        }

        for(auto inclusion : this->_inclusions){
            inclusion->eshelby_tensor()->reinit();
        }

        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();

        //std::cout<<"Start "<<std::endl;
        while (true) {
            set_local_history();
            set_initial_values();

            //update_strain();
            for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
                dFi[i] = tmech::dcontract(this->_strain_concentration_tensors[i], this->_dstrain);
                if(dynamic_cast<incremental_solid_material<_T, _Dim>*>(_finite_strain_materials[i])){
                    dynamic_cast<incremental_solid_material<_T, _Dim>*>(_finite_strain_materials[i])->dstrain_tensor() = dFi[i];
                }else{
                    _finite_strain_materials[i]->deformation_tensor() += dFi[i];
                }
            }

            for(auto solid : _solid_materials){
                solid->update();
            }

            for(auto inclusion : this->_inclusions){
                inclusion->eshelby_tensor()->reinit();
            }

            update_volume_fraction();

            //update strain concentratoin tensors
            mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();

            //update_R(R);
            for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
                R[i] = dFi[i] - tmech::dcontract(this->_strain_concentration_tensors[i], this->_dstrain);
            }

            const value_type norm = R_norm(R);
            //            if(norm != 0){
            //                std::cout<<"Norm: "<<norm<<std::endl;
            //            }
            //std::cout<<"Norm: "<<norm<<std::endl;
            if(norm<=tol){break;}

            if(iter == max_iter){
                std::cout<<"Norm: "<<norm<<std::endl;
                throw std::runtime_error("finite_strain_mean_field_composite::update() max number of iteration reached");
            }
            ++iter;
        }

        for(auto solid : _solid_materials){
            solid->update();
        }

        update_macro_stress();
        update_macro_tangent();

        get_local_history();
        copy_illustrative_to_history();
    }

    inline virtual void update_stress() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_stress(): is not initialized");
        }

        update();
        //        set_local_history_all();

        //        update_strain();

        //        for(auto solid : _solid_materials){
        //            solid->update_stress();
        //        }
        //        update_macro_stress();

        //        get_local_history_all();
    }

    inline virtual void update_tangent() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_tangent(): is not initialized");
        }

        update();
        //        update_strain();
        //        for(auto solid : _solid_materials){
        //            solid->update_tangent();
        //        }
        //        update_macro_tangent();
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override {
        return this;
    }


private:
    constexpr inline auto update_R(std::vector<tensor2> & __R)const{
        //        __R.back() = _matrix_material.get()->dstrain() - tmech::dcontract(_strain_ct.back(), this->_dstrain);
        //        for(size_type i{0}; i<_inclusion_material.size(); ++i){
        //            __R[i] = _inclusion_material[i].get()->dstrain() - tmech::dcontract(_strain_ct[i], this->_dstrain);
        //        }
    }

    constexpr inline auto R_norm(std::vector<tensor2> const& __R)const{
        value_type norm{0};
        for(const auto& R : __R){
            const auto temp = tmech::norm(R);
            norm += temp*temp;
        }
        return std::sqrt(norm);
    }

    constexpr inline auto update_volume_fraction(){
        const auto detJ_macro{tmech::det(this->deformation_tensor() + this->_dstrain)};
        value_type sum{0};
        for(auto [c, solid, finite, A] : _reference_formulations){
            const auto F{finite->deformation_tensor()};
            c = tmech::det(F)*c/detJ_macro;
            sum += c;
            //std::cout<<c<<std::endl;
        }

        for(auto [c, solid, finite, A] : _mixed_formulations){
            const auto F{finite->deformation_tensor()};
            c = tmech::det(F)*c/detJ_macro;
            sum += c;
        }

        //std::cout<<1/sum<<" "<<detJ_macro<<std::endl;

        for(auto [c, solid, finite, A] : _reference_formulations){
            c /= sum;
        }
        for(auto [c, solid, finite, A] : _mixed_formulations){
            c /= sum;
        }
    }

    constexpr inline auto copy_illustrative_to_history(){
        const auto start{history_base::_history.size() - _init_values.size()};
        size_type iter{0};
        history_base::_history[start + iter++] = this->_cm;
        for(auto inclusion : this->_inclusions){
            history_base::_history[start + iter++] = inclusion->volume_fraction();
        }
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    history_base::_history[start + iter++] = S->parameter()[i];
                }
            }
        }
    }

    constexpr inline auto get_initial_values(){
        size_type iter{0};
        //volume fractions
        _init_values[iter++] = this->_cm;
        for(auto inclusion : this->_inclusions){
            _init_values[iter++] = inclusion->volume_fraction();
        }
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    _init_values[iter++] = S->parameter()[i];
                }
            }
        }
    }

    constexpr inline auto set_initial_values(){
        size_type iter{0};
        //volume fractions
        this->_cm = _init_values[iter++];
        for(auto inclusion : this->_inclusions){
            inclusion->volume_fraction() = _init_values[iter++];
        }
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    S->parameter()[i] = _init_values[iter++];
                }
            }
        }
    }

    constexpr inline auto get_local_history(){
        history_base::get_local_history();
        size_type iter{0};
        //set deformation tensor
        for(auto material : _finite_strain_materials){
            adaptor(&history_base::_history[_local_history_size + iter]) = material->deformation_tensor();
            iter += _Dim*_Dim;
        }
    }

    constexpr inline auto set_local_history(){
        history_base::set_local_history();
        size_type iter{0};
        //set deformation tensor
        for(auto material : _finite_strain_materials){
            material->deformation_tensor() = adaptor(&history_base::_history[_local_history_size + iter]);
            iter += _Dim*_Dim;
        }
    }

    //    constexpr inline auto get_local_history_all(){
    //        history_base::get_local_history();
    //        size_type iter{0};
    //        //set deformation tensor
    //        for(auto material : _finite_strain_materials){
    //            adaptor(&history_base::_history[_local_history_size + iter]) = material->deformation_tensor();
    //            iter += (_Dim == 2 ? 3 : 6);
    //        }
    //        //volume fractions
    //        this->_history[_local_history_size + iter++] = this->_cm;
    //        for(auto inclusion : this->_inclusions){
    //            history_base::_history[_local_history_size + iter++] = inclusion->volume_fraction();
    //        }
    //        //geometry eshelby tensors
    //        for(auto inclusion : this->_inclusions){
    //            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
    //            if(S){
    //                for(size_type i{0}; i<4; ++i){
    //                    history_base::_history[_local_history_size + iter++] = S->parameter()[i];
    //                }
    //            }
    //        }
    //        //rotation eshelby tensors
    //        for(auto inclusion : this->_inclusions){
    //            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
    //            if(S){
    //                //adaptor(&history_base::_history[_local_history_size + iter]) = S->rotation_tensor();
    //                iter += (_Dim == 2 ? 3 : 6);
    //            }
    //        }
    //    }


    //    constexpr inline auto set_local_history_all(){
    //        history_base::set_local_history();
    //        size_type iter{0};
    //        //set deformation tensor
    //        for(auto material : _finite_strain_materials){
    //            material->deformation_tensor() = adaptor(&history_base::_history[_local_history_size + iter]);
    //            iter += (_Dim == 2 ? 3 : 6);
    //        }
    //        //volume fractions
    //        this->_cm = this->_history[_local_history_size + iter++];
    //        for(auto inclusion : this->_inclusions){
    //            inclusion->volume_fraction() = history_base::_history[_local_history_size + iter++];
    //        }
    //        //geometry eshelby tensors
    //        for(auto inclusion : this->_inclusions){
    //            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
    //            if(S){
    //                //std::cout<<"History ";
    //                for(size_type i{0}; i<4; ++i){
    //                    S->parameter()[i] = history_base::_history[_local_history_size + iter++];
    //                    //std::cout<<S->parameter()[i]<<" ";
    //                }
    //                //std::cout<<std::endl;
    //            }
    //        }
    //        //rotation eshelby tensors
    //        for(auto inclusion : this->_inclusions){
    //            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
    //            if(S){
    //                //S->rotation_tensor() = adaptor(&history_base::_history[_local_history_size + iter]);
    //                iter += (_Dim == 2 ? 3 : 6);
    //            }
    //        }
    //    }

    constexpr inline auto update_macro_stress(){
        this->_stress.fill(0);
        update_macro_stress_reference();
        update_macro_stress_mixed();
    }

    constexpr inline auto update_macro_stress_reference(){
        for(auto [c, solid, finite, A] : _reference_formulations){
            const auto& F{finite->deformation_tensor()};
            this->_stress += c*F*solid->stress_tensor();
        }
    }

    constexpr inline auto update_macro_stress_mixed(){
        for(auto [c, solid, finit, A] : _mixed_formulations){
            this->_stress += c*solid->stress_tensor();
        }
    }


    constexpr inline auto update_macro_tangent(){
        this->_C.fill(0);
        update_macro_tangent_reference();
        update_macro_tangent_mixed();
    }

    constexpr inline auto update_macro_tangent_reference(){
        const tmech::eye<value_type, _Dim, 2> I;
        for(auto [c, solid, finite, A] : _reference_formulations){
            const auto& F{finite->deformation_tensor()};
            this->_C += c*tmech::dcontract(tmech::dcontract(tmech::dcontract(tmech::otimesu(F, I), solid->tangent_tensor()), tmech::otimesu(tmech::trans(F), I))
                                           + tmech::otimesu(I, solid->stress_tensor()), A);
        }
    }

    constexpr inline auto update_macro_tangent_mixed(){
        for(auto [c, solid, finit, A] : _mixed_formulations){
            this->_C += c*tmech::dcontract(solid->tangent_tensor(), A);
        }
    }

    std::vector<finite_strain_solid_material_base<_T, _Dim, _Container>*> _finite_strain_materials;
    std::vector<solid_material_base<_T, _Dim, _Container>*> _solid_materials;
    std::vector<mixed_formulation_wrapper> _mixed_formulations;
    std::vector<reference_formulation_wrapper> _reference_formulations;
    size_type _local_history_size;
    std::vector<value_type> _init_values;
};








template <typename _T, std::size_t _Dim, typename _Container>
class finite_strain_mean_field_composite_explicit :
        public solid_material_base<_T, _Dim, _Container>,
        public mean_field_composite_solid_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>,
        public incremental_solid_material<_T, _Dim>,
        public composite_history_base<_T, _Dim, _Container>
{
    //Mixed <vol_frac, solid_material, strain_ct>
    struct mixed_formulation_wrapper
    {
        _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        finite_strain_solid_material_base<_T, _Dim, _Container>* _finite_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
    //Reference <vol_frac, solid_material, finite_material, strain_ct>
    struct reference_formulation_wrapper
    {
        _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        finite_strain_solid_material_base<_T, _Dim, _Container>* _finite_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
    using adaptor = tmech::adaptor<_T,_Dim, 2, tmech::full<_Dim>>;
    using adaptor4 = tmech::adaptor<_T,_Dim, 4, tmech::full<_Dim>>;
    using tensor2 = tmech::tensor<_T, _Dim, 2>;
public:
    using value_type = _T;
    using size_type  = std::size_t;
    using history_base = composite_history_base<_T, _Dim, _Container>;

    finite_strain_mean_field_composite_explicit():
        solid_material_base<_T, _Dim, _Container>(),
        mean_field_composite_solid_base<_T, _Dim, _Container>(),
        finite_strain_solid_material_base<_T, _Dim, _Container>(this),
        incremental_solid_material<_T, _Dim>(),
        composite_history_base<_T, _Dim, _Container>(this),
        _finite_strain_materials(),
        _solid_materials(),
        _mixed_formulations(),
        _reference_formulations(),
        _local_history_size(0)//,
      //_macro_vol_frac(1)
    {}

    virtual ~finite_strain_mean_field_composite_explicit() {}

    inline virtual void init() override {
        if(!solid_material_base<_T, _Dim, _Container>::_is_init || !mean_field_composite_solid_base<_T, _Dim, _Container>::_is_init){
            this->_matrix_material->init();
            for(auto inclusion : this->_inclusions){
                inclusion->init();
            }

            finite_strain_solid_material_base<_T, _Dim, _Container>::_type = FINITE_STRAIN_FORMULATION::MIXED;

            mean_field_composite_solid_base<_T, _Dim, _Container>::init();

            //setup finite strain material pointer
            _finite_strain_materials.resize(this->_inclusions.size()+1);
            _finite_strain_materials.back() = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _finite_strain_materials[i] = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }

            //setup solid material pointer
            _solid_materials.resize(this->_inclusions.size()+1);
            _solid_materials.back() = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _solid_materials[i] = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }

            //setup mixed formulation wrappers
            //and init F_i = I
            for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
                _finite_strain_materials[i]->deformation_tensor() = tmech::eye<value_type,_Dim,2>();
                auto& vol_frac{(i == _finite_strain_materials.size() - 1 ? this->_cm : this->_inclusions[i]->volume_fraction())};
                switch (_finite_strain_materials[i]->formulation_type()) {
                case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
                    _reference_formulations.push_back({vol_frac, _solid_materials[i], _finite_strain_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                case FINITE_STRAIN_FORMULATION::MIXED:
                    _mixed_formulations.push_back({vol_frac, _solid_materials[i], _finite_strain_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                }
            }

            //init history
            //get local history size
            history_base::init();
            _local_history_size = history_base::size();
            auto his_size{history_base::size()};
            //entries of each deformation tensor
            his_size += _Dim*_Dim*_solid_materials.size();
            //entries of each strain concentration tensor
            his_size += _Dim*_Dim*_Dim*_Dim*_solid_materials.size();
            //volume fraction for each material phase
            his_size += _solid_materials.size();

            //initial values
            //vol fraction at reference solution
            auto init_val{_solid_materials.size()};
            //for each finite_strain_eshelby_tensor
            //-->  geometry at tn
            for(auto inclusion : this->_inclusions){
                auto S = inclusion->eshelby_tensor();
                if(dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(S)){
                    init_val += 3;//geometry + nu
                }
            }

            //get init volume fraction, geometry
            _init_values.resize(init_val+1);
            get_initial_values();

            //resize history
            history_base::resize(his_size+init_val);//init_val only for illustration + 1 macro volume fraction

            get_local_history();

            update_macro_tangent();

            solid_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    inline virtual void reinit() override {
        this->_matrix_material->reinit();
        for(auto inclusion : this->_inclusions){
            inclusion->reinit();
        }
        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();
        update_macro_tangent();
    }

    inline virtual void update_strain()override{
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_strain(): class is not initialized");
        }

        //const auto& F_macro{this->deformation_tensor() + this->_dstrain};
        for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
            auto incremental{dynamic_cast<incremental_solid_material<_T,_Dim>*>(_finite_strain_materials[i])};
            if(incremental){
                incremental->dstrain_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], this->_dstrain);
                _finite_strain_materials[i]->deformation_tensor_n() = _finite_strain_materials[i]->deformation_tensor();
                _finite_strain_materials[i]->deformation_tensor() += incremental->dstrain_tensor();
            }else{
                _finite_strain_materials[i]->deformation_tensor_n() = _finite_strain_materials[i]->deformation_tensor();
                _finite_strain_materials[i]->deformation_tensor() += tmech::dcontract(this->_strain_concentration_tensors[i], this->_dstrain);
            }
        }
    }

    inline void update_material(){
        //update material
        for(auto solid : _solid_materials){
            auto finite_strain = dynamic_cast<finite_strain_mean_field_composite_explicit<_T,_Dim,_Container>*>(solid);
            if(finite_strain){
                finite_strain->set_local_history();
                finite_strain->update_strain();
                finite_strain->update_material();
                finite_strain->update_macroscopic_properties();
            }else{
                solid->update();
            }
        }
    }

    inline void update_strain_concentration_tensor(){
        for(auto solid : _solid_materials){
            auto finite_strain = dynamic_cast<finite_strain_mean_field_composite_explicit<_T,_Dim,_Container>*>(solid);
            if(finite_strain){
                finite_strain->update_strain_concentration_tensor();
            }
        }

        //update eshely tensor with F^{n+1} for next step
        for(auto inclusion : this->_inclusions){
            inclusion->eshelby_tensor()->reinit();
        }

        //update volume fraction using F^{n+1}
        //update_volume_fraction();

        //update strain concentration tensors with quantities from last step
        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();
    }

    inline void update_macroscopic_properties(){
        this->_stress.fill(0);
        this->_C.fill(0);
        for(std::size_t i{0};i<_finite_strain_materials.size();++i){
            const auto c{this->volume_fraction(i)};
            //std::cout<<c<<std::endl;
            //std::cout<<c<<std::endl;
            this->_stress += c*_finite_strain_materials[i]->mixed_stress_tensor();
            this->_C += c*tmech::dcontract(_finite_strain_materials[i]->mixed_tangent_tensor(),this->_strain_concentration_tensors[i]);
            //std::cout<<"material "<<i<<"\n"<<_finite_strain_materials[i]->mixed_tangent_tensor()<<std::endl;
        }
    }

    inline virtual void update() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update(): is not initialized");
        }

        //set history and deformation tensor t^n
        set_local_history();

        update_strain();

        update_material();

        //initial values of S and volume fractions
        set_initial_values();

        update_macroscopic_properties();

        update_strain_concentration_tensor();

        get_local_history();
        //copy_illustrative_to_history();
    }

    inline virtual void update_stress() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_stress(): is not initialized");
        }

        set_local_history();
        set_initial_values();

        //update deformation tensor F^{n+1}, this is not done at the micro level
        for(std::size_t i{0}; i<_solid_materials.size(); ++i){
            auto incremental{dynamic_cast<incremental_solid_material<_T,_Dim>*>(_solid_materials[i])};
            if(incremental){
                incremental->dstrain_tensor().fill(0);
            }
        }

        //call material with deformation tensor F^n
        for(auto solid : _solid_materials){
            //dynamic_cast<incremental_solid_material<_T,_Dim>*>(solid)->dstrain_tensor().fill(0);
            solid->update();
        }

        //update eshely tensor with F^n
        for(auto inclusion : this->_inclusions){
            //inclusion->eshelby_tensor()->reinit();
        }

        //update volume fraction using F^n
        //update_volume_fraction();

        //update strain concentration tensors with quantities from last step
        //mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();

        //set history and deformation tensor last step
        set_local_history();

        update_strain();

        //update material
        for(auto solid : _solid_materials){
            solid->update_stress();
        }

        update_macro_stress();

        get_local_history();
        copy_illustrative_to_history();
    }

    inline virtual void update_tangent() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_tangent(): is not initialized");
        }

        update();
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override {
        return this;
    }

    constexpr inline auto update_volume_fraction(){
        //const auto detJ_macro{tmech::det(this->deformation_tensor())};
        value_type sum{0};
        for(auto [c, solid, finite, A] : _reference_formulations){
            const auto F{finite->deformation_tensor()};//is still the old one
            c = tmech::det(F)*c;///detJ_macro;
            sum += c;
        }

        for(auto [c, solid, finite, A] : _mixed_formulations){
            const auto F{finite->deformation_tensor()};
            c = tmech::det(F)*c;///detJ_macro;
            sum += c;
        }

        for(auto [c, solid, finite, A] : _reference_formulations){
            c /= sum;

        }
        for(auto [c, solid, finite, A] : _mixed_formulations){
            c /= sum;
        }

        //std::cout<<this->volume_fraction(0)<<" "<<this->volume_fraction(1)<<" "<<this->volume_fraction(0)+this->volume_fraction(1)<<std::endl;
    }

    constexpr inline auto get_initial_values(){
        size_type iter{0};
        //volume fractions
        //_init_values[iter++] = this->_cm;
        for(auto inclusion : this->_inclusions){
            //_init_values[iter++] = inclusion->volume_fraction();
        }
        //macro volume fraction
        //_init_values[iter++] = this->_parameter[0];
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    _init_values[iter++] = S->parameter()[i];
                }
            }
        }
    }

    constexpr inline void set_initial_values(){
        size_type iter{0};
        //volume fractions
        //this->_cm = _init_values[iter++];
        for(auto inclusion : this->_inclusions){
            //inclusion->volume_fraction() = _init_values[iter++];
        }
        //macro volume fraction
        //this->_parameter[0] = _init_values[iter++];
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    S->parameter()[i] = _init_values[iter++];
                }
            }
        }
    }

    constexpr inline void get_local_history(){

        for(auto solid : _solid_materials){
            auto finite_strain = dynamic_cast<finite_strain_mean_field_composite_explicit<_T,_Dim,_Container>*>(solid);
            if(finite_strain){
                finite_strain->get_local_history();
            }
        }

        history_base::get_local_history();
        size_type iter{0};
        //set deformation tensor
        for(std::size_t i{0}; i<_solid_materials.size(); ++i){
            adaptor(&history_base::_history[_local_history_size + iter]) = _finite_strain_materials[i]->deformation_tensor();
            iter += _Dim*_Dim;
            //            adaptor(&history_base::_history[_local_history_size + iter]) = _solid_materials[i]->stress_tensor();
            //            iter += _Dim*_Dim;
        }
        for(auto const& A : this->_strain_concentration_tensors){
            adaptor4(&history_base::_history[_local_history_size + iter]) = A;
            iter += _Dim*_Dim*_Dim*_Dim;
        }

        //        for(std::size_t i{0}; i<_solid_materials.size(); ++i){
        //            history_base::_history[_local_history_size + iter++] = this->volume_fraction(i);
        //        }
    }

    constexpr inline auto set_local_history(){
        history_base::set_local_history();
        size_type iter{0};
        //set deformation tensor
        for(std::size_t i{0}; i<_solid_materials.size(); ++i){
            _finite_strain_materials[i]->deformation_tensor() = adaptor(&history_base::_history[_local_history_size + iter]);
            iter += _Dim*_Dim;
            //            _solid_materials[i]->stress_tensor() = adaptor(&history_base::_history[_local_history_size + iter]);
            //            iter += _Dim*_Dim;
        }

        for(auto& A : this->_strain_concentration_tensors){
            A = adaptor4(&history_base::_history[_local_history_size + iter]);
            iter += _Dim*_Dim*_Dim*_Dim;
        }

        //        for(std::size_t i{0}; i<_solid_materials.size(); ++i){
        //            this->volume_fraction(i) = history_base::_history[_local_history_size + iter++];
        //        }
    }

    //    constexpr inline value_type& macro_volume_fraction(){
    //        return _macro_vol_frac;
    //    }


private:
    constexpr inline auto copy_illustrative_to_history(){
        const auto start{history_base::_history.size() - _init_values.size()};
        size_type iter{0};
        history_base::_history[start + iter++] = this->_cm;
        for(auto inclusion : this->_inclusions){
            history_base::_history[start + iter++] = inclusion->volume_fraction();
        }
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    history_base::_history[start + iter++] = S->parameter()[i];
                }
            }
        }
    }

    constexpr inline auto update_macro_stress(){
        this->_stress.fill(0);
        update_macro_stress_reference();
        update_macro_stress_mixed();
    }

    constexpr inline auto update_macro_stress_reference(){
        for(auto [c, solid, finite, A] : _reference_formulations){
            //const auto& F{finite->deformation_tensor()};
            //this->_stress += c*F*solid->stress_tensor();//
            this->_stress += c*finite->mixed_stress_tensor();//
        }
    }

    constexpr inline auto update_macro_stress_mixed(){
        for(auto [c, solid, finite, A] : _mixed_formulations){
            this->_stress += c*finite->mixed_stress_tensor();
        }
    }


    constexpr inline auto update_macro_tangent(){
        this->_C.fill(0);
        update_macro_tangent_reference();
        update_macro_tangent_mixed();
    }

    constexpr inline auto update_macro_tangent_reference(){
        const tmech::eye<value_type, _Dim, 2> I;
        for(auto [c, solid, finite, A] : _reference_formulations){
            //const auto& F{finite->deformation_tensor()};
            this->_C += c*tmech::dcontract(finite->mixed_tangent_tensor(),A);//
            //            this->_C += tmech::dcontract(tmech::dcontract(tmech::dcontract(tmech::otimesu(F, I), solid->tangent_tensor()), tmech::otimesu(tmech::trans(F), I))
            //                            + tmech::otimesu(I, solid->stress_tensor()), A);
            //this->_C += tmech::dcontract(tmech::dcontract(tmech::otimesu(F, I), solid->tangent_tensor()), tmech::otimesu(tmech::trans(F), I))
            //                + tmech::otimesu(I, solid->stress_tensor());
        }
    }

    constexpr inline auto update_macro_tangent_mixed(){
        for(auto [c, solid, finite, A] : _mixed_formulations){
            this->_C += c*tmech::dcontract(finite->mixed_tangent_tensor(), A);
        }
    }

    std::vector<finite_strain_solid_material_base<_T, _Dim, _Container>*> _finite_strain_materials;
    std::vector<solid_material_base<_T, _Dim, _Container>*> _solid_materials;
    std::vector<mixed_formulation_wrapper> _mixed_formulations;
    std::vector<reference_formulation_wrapper> _reference_formulations;
    size_type _local_history_size;
    std::vector<value_type> _init_values;
    //value_type _macro_vol_frac;
};




template <typename _T, std::size_t _Dim, typename _Container>
class finite_strain_mean_field_composite_explicit_incremental :
        public mean_field_composite_solid_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>,
        public incremental_solid_material<_T, _Dim>,
        public composite_history_base<_T, _Dim, _Container>,
        public solid_material_base<_T, _Dim, _Container>
{
    //Mixed <vol_frac, solid_material, strain_ct>
    struct mixed_formulation_wrapper
    {
        _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        finite_strain_solid_material_base<_T, _Dim, _Container>* _finite_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
    //Reference <vol_frac, solid_material, finite_material, strain_ct>
    struct reference_formulation_wrapper
    {
        _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        finite_strain_solid_material_base<_T, _Dim, _Container>* _finite_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
    using adaptor = tmech::adaptor<_T,_Dim, 2, tmech::voigt<_Dim>>;
    using tensor2 = tmech::tensor<_T, _Dim, 2>;
public:
    using value_type = _T;
    using size_type  = std::size_t;
    using history_base = composite_history_base<_T, _Dim, _Container>;

    finite_strain_mean_field_composite_explicit_incremental():
        mean_field_composite_solid_base<_T, _Dim, _Container>(),
        finite_strain_solid_material_base<_T, _Dim, _Container>(this),
        incremental_solid_material<_T, _Dim>(),
        composite_history_base<_T, _Dim, _Container>(this),
        _finite_strain_materials(),
        _solid_materials(),
        _mixed_formulations(),
        _reference_formulations(),
        _local_history_size(0)
    {}

    virtual ~finite_strain_mean_field_composite_explicit_incremental() {}

    inline virtual void init() override {
        if(!solid_material_base<_T, _Dim, _Container>::_is_init || !mean_field_composite_solid_base<_T, _Dim, _Container>::_is_init){
            this->_matrix_material->init();
            for(auto inclusion : this->_inclusions){
                inclusion->init();
            }

            finite_strain_solid_material_base<_T, _Dim, _Container>::_type = FINITE_STRAIN_FORMULATION::MIXED;

            mean_field_composite_solid_base<_T, _Dim, _Container>::init();
            update_macro_tangent();

            //setup finite strain material pointer
            _finite_strain_materials.resize(this->_inclusions.size()+1);
            _finite_strain_materials.back() = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _finite_strain_materials[i] = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }

            //setup solid material pointer
            _solid_materials.resize(this->_inclusions.size()+1);
            _solid_materials.back() = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _solid_materials[i] = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }

            //setup mixed formulation wrappers
            //and init F_i = I
            for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
                _finite_strain_materials[i]->deformation_tensor() = tmech::eye<value_type,_Dim,2>();
                auto& vol_frac{(i == _finite_strain_materials.size() - 1 ? this->_cm : this->_inclusions[i]->volume_fraction())};
                switch (_finite_strain_materials[i]->formulation_type()) {
                case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
                    _reference_formulations.push_back({vol_frac, _solid_materials[i], _finite_strain_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                case FINITE_STRAIN_FORMULATION::MIXED:
                    _mixed_formulations.push_back({vol_frac, _solid_materials[i], _finite_strain_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                }
            }

            //init history
            //get local history size
            history_base::init();
            _local_history_size = history_base::size();
            auto his_size{history_base::size()};
            //entries of each deformation tensor
            his_size += (_Dim == 2 ? 3 : 6)*_solid_materials.size()*2;

            //initial values
            //vol fraction at reference solution
            auto init_val{_solid_materials.size()};
            //for each finite_strain_eshelby_tensor
            //-->  geometry at tn
            for(auto inclusion : this->_inclusions){
                auto S = inclusion->eshelby_tensor();
                if(dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(S)){
                    init_val += 3;//geometry + nu
                }
            }

            //get init volume fraction, geometry
            _init_values.resize(init_val);
            get_initial_values();

            //resize history
            history_base::resize(his_size+init_val);//init_val only for illustration

            get_local_history();

            solid_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    inline virtual void reinit() override {
        this->_matrix_material->reinit();
        for(auto inclusion : this->_inclusions){
            inclusion->reinit();
        }
        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();
        update_macro_tangent();
    }

    inline virtual void update_strain()override{
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_strain(): class is not initialized");
        }

        //const auto& F_macro{this->deformation_tensor() + this->_dstrain};
        for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
            auto incremental{dynamic_cast<incremental_solid_material<_T,_Dim>*>(_finite_strain_materials[i])};
            if(incremental){
                incremental->dstrain_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], this->_dstrain);
                _finite_strain_materials[i]->deformation_tensor_n() = _finite_strain_materials[i]->deformation_tensor();
                _finite_strain_materials[i]->deformation_tensor() += incremental->dstrain_tensor();
            }else{
                _finite_strain_materials[i]->deformation_tensor_n() = _finite_strain_materials[i]->deformation_tensor();
                _finite_strain_materials[i]->deformation_tensor() += tmech::dcontract(this->_strain_concentration_tensors[i], this->_dstrain);
            }
        }
    }

    inline virtual void update() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update(): is not initialized");
        }

        //std::cout<<"Start update"<<std::endl;
        //set history and deformation tensor last step
        set_local_history();
        set_initial_values();

        for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
            auto incremental{dynamic_cast<incremental_solid_material<_T,_Dim>*>(_finite_strain_materials[i])};
            if(incremental){
                incremental->dstrain_tensor().fill(0);
            }
        }

        //call material with deformation tensor F^n
        for(auto solid : _solid_materials){
            solid->update();
            //std::cout<<solid->stress_tensor()<<std::endl;
        }

        //update eshely tensor with F^n
        for(auto inclusion : this->_inclusions){
            inclusion->eshelby_tensor()->reinit();
        }

        //update volume fraction using F^n
        update_volume_fraction();

        //update strain concentration tensors with quantities from last step
        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();

        update_strain();


        //update material
        for(auto solid : _solid_materials){
            solid->update();
        }

        //this->_stress.fill(0);
        this->_C.fill(0);
        for(std::size_t i{0};i<_finite_strain_materials.size();++i){
            const auto c{this->volume_fraction(i)};
            this->_stress += c*_finite_strain_materials[i]->mixed_stress_tensor();
            this->_C += c*tmech::dcontract(_finite_strain_materials[i]->mixed_tangent_tensor(),this->_strain_concentration_tensors[i]);
        }

        //update_macro_stress();
        //update_macro_tangent();

        get_local_history();
        copy_illustrative_to_history();
    }

    inline virtual void update_stress() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_stress(): is not initialized");
        }

        set_local_history();
        set_initial_values();

        //update deformation tensor F^{n+1}, this is not done at the micro level
        for(std::size_t i{0}; i<_solid_materials.size(); ++i){
            auto incremental{dynamic_cast<incremental_solid_material<_T,_Dim>*>(_solid_materials[i])};
            if(incremental){
                incremental->dstrain_tensor().fill(0);
            }
        }

        //call material with deformation tensor F^n
        for(auto solid : _solid_materials){
            //dynamic_cast<incremental_solid_material<_T,_Dim>*>(solid)->dstrain_tensor().fill(0);
            solid->update();
        }

        //update eshely tensor with F^n
        for(auto inclusion : this->_inclusions){
            //inclusion->eshelby_tensor()->reinit();
        }

        //update volume fraction using F^n
        //update_volume_fraction();

        //update strain concentration tensors with quantities from last step
        //mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();

        //set history and deformation tensor last step
        set_local_history();

        update_strain();

        //update material
        for(auto solid : _solid_materials){
            solid->update_stress();
        }

        update_macro_stress();

        get_local_history();
        copy_illustrative_to_history();
    }

    inline virtual void update_tangent() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_tangent(): is not initialized");
        }

        update();
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override {
        return this;
    }


private:

    constexpr inline auto update_volume_fraction(){
        const auto detJ_macro{tmech::det(this->deformation_tensor_n())};
        value_type sum{0};
        for(auto [c, solid, finite, A] : _reference_formulations){
            const auto F{finite->deformation_tensor()};//is still the old one
            c = tmech::det(F)*c/detJ_macro;
            sum += c;
        }

        for(auto [c, solid, finite, A] : _mixed_formulations){
            const auto F{finite->deformation_tensor()};
            c = tmech::det(F)*c/detJ_macro;
            sum += c;
        }

        for(auto [c, solid, finite, A] : _reference_formulations){
            c /= sum;
        }
        for(auto [c, solid, finite, A] : _mixed_formulations){
            c /= sum;
        }
    }

    constexpr inline auto copy_illustrative_to_history(){
        const auto start{history_base::_history.size() - _init_values.size()};
        size_type iter{0};
        history_base::_history[start + iter++] = this->_cm;
        for(auto inclusion : this->_inclusions){
            history_base::_history[start + iter++] = inclusion->volume_fraction();
        }
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    history_base::_history[start + iter++] = S->parameter()[i];
                }
            }
        }
    }

    constexpr inline auto get_initial_values(){
        size_type iter{0};
        //volume fractions
        _init_values[iter++] = this->_cm;
        for(auto inclusion : this->_inclusions){
            _init_values[iter++] = inclusion->volume_fraction();
        }
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    _init_values[iter++] = S->parameter()[i];
                }
            }
        }
    }

    constexpr inline auto set_initial_values(){
        size_type iter{0};
        //volume fractions
        this->_cm = _init_values[iter++];
        for(auto inclusion : this->_inclusions){
            inclusion->volume_fraction() = _init_values[iter++];
        }
        //geometry eshelby tensors
        for(auto inclusion : this->_inclusions){
            auto S = dynamic_cast<eshelby_tensor_finite_strain_base<_T,_Dim,_Container>*>(inclusion->eshelby_tensor());
            if(S){
                for(size_type i{0}; i<3; ++i){
                    S->parameter()[i] = _init_values[iter++];
                }
            }
        }
    }

    constexpr inline auto get_local_history(){
        history_base::get_local_history();
        size_type iter{0};
        //set deformation tensor
        for(auto material : _finite_strain_materials){
            adaptor(&history_base::_history[_local_history_size + iter]) = material->deformation_tensor();
            iter += (_Dim == 2 ? 3 : 6);
        }
        for(auto material : _solid_materials){
            adaptor(&history_base::_history[_local_history_size + iter]) = material->stress_tensor();
            iter += (_Dim == 2 ? 3 : 6);
        }
    }

    constexpr inline auto set_local_history(){
        history_base::set_local_history();
        size_type iter{0};
        //set deformation tensor
        for(auto material : _finite_strain_materials){
            material->deformation_tensor() = adaptor(&history_base::_history[_local_history_size + iter]);
            iter += (_Dim == 2 ? 3 : 6);
        }
        for(auto material : _solid_materials){
            material->stress_tensor() = adaptor(&history_base::_history[_local_history_size + iter]);
            iter += (_Dim == 2 ? 3 : 6);
        }
    }

    constexpr inline auto update_macro_stress(){
        this->_stress.fill(0);
        update_macro_stress_reference();
        update_macro_stress_mixed();
    }

    constexpr inline auto update_macro_stress_reference(){
        for(auto [c, solid, finite, A] : _reference_formulations){
            //const auto& F{finite->deformation_tensor()};
            //this->_stress += c*F*solid->stress_tensor();//
            this->_stress += c*finite->mixed_stress_tensor();//
        }
    }

    constexpr inline auto update_macro_stress_mixed(){
        for(auto [c, solid, finite, A] : _mixed_formulations){
            this->_stress += c*finite->mixed_stress_tensor();
        }
    }


    constexpr inline auto update_macro_tangent(){
        this->_C.fill(0);
        update_macro_tangent_reference();
        update_macro_tangent_mixed();
    }

    constexpr inline auto update_macro_tangent_reference(){
        const tmech::eye<value_type, _Dim, 2> I;
        for(auto [c, solid, finite, A] : _reference_formulations){
            //const auto& F{finite->deformation_tensor()};
            this->_C += c*tmech::dcontract(finite->mixed_tangent_tensor(),A);//
            //            this->_C += tmech::dcontract(tmech::dcontract(tmech::dcontract(tmech::otimesu(F, I), solid->tangent_tensor()), tmech::otimesu(tmech::trans(F), I))
            //                            + tmech::otimesu(I, solid->stress_tensor()), A);
            //this->_C += tmech::dcontract(tmech::dcontract(tmech::otimesu(F, I), solid->tangent_tensor()), tmech::otimesu(tmech::trans(F), I))
            //                + tmech::otimesu(I, solid->stress_tensor());
        }
    }

    constexpr inline auto update_macro_tangent_mixed(){
        for(auto [c, solid, finite, A] : _mixed_formulations){
            this->_C += c*tmech::dcontract(finite->mixed_tangent_tensor(), A);
        }
    }

    std::vector<finite_strain_solid_material_base<_T, _Dim, _Container>*> _finite_strain_materials;
    std::vector<solid_material_base<_T, _Dim, _Container>*> _solid_materials;
    std::vector<mixed_formulation_wrapper> _mixed_formulations;
    std::vector<reference_formulation_wrapper> _reference_formulations;
    size_type _local_history_size;
    std::vector<value_type> _init_values;
};




template <typename _T, std::size_t _Dim, typename _Container>
class finite_strain_mean_field_composite_elastic :
        public mean_field_composite_solid_base<_T, _Dim, _Container>,
        public finite_strain_solid_material_base<_T, _Dim, _Container>,
        public solid_material_base<_T, _Dim, _Container>
        //public incremental_solid_material<_T, _Dim>//,
        //public composite_history_base<_T, _Dim, _Container>
{
    //Mixed <vol_frac, solid_material, strain_ct>
    struct mixed_formulation_wrapper
    {
        const _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
    //Reference <vol_frac, solid_material, finite_material, strain_ct>
    struct reference_formulation_wrapper
    {
        const _T& vol_frac;
        solid_material_base<_T, _Dim, _Container>* _solid_material;
        finite_strain_solid_material_base<_T, _Dim, _Container>* _finite_material;
        tmech::tensor<_T, _Dim, 4> const& _A;
    };
public:
    using value_type = _T;
    using size_type  = std::size_t;
    using history_base = composite_history_base<_T, _Dim, _Container>;

    finite_strain_mean_field_composite_elastic():
        mean_field_composite_solid_base<_T, _Dim, _Container>(),
        finite_strain_solid_material_base<_T, _Dim, _Container>(this)//,
      //incremental_solid_material<_T, _Dim>()//,
      //composite_history_base<_T, _Dim, _Container>(this)
    {}

    virtual ~finite_strain_mean_field_composite_elastic() {}

    inline virtual void init() override {
        if(!solid_material_base<_T, _Dim, _Container>::_is_init || !mean_field_composite_solid_base<_T, _Dim, _Container>::_is_init){
            this->_matrix_material->init();
            for(auto inclusion : this->_inclusions){
                inclusion->init();
            }

            finite_strain_solid_material_base<_T, _Dim, _Container>::_type = FINITE_STRAIN_FORMULATION::MIXED;

            mean_field_composite_solid_base<_T, _Dim, _Container>::init();
            update_macro_tangent();

            //setup finite strain material pointer
            _finite_strain_materials.resize(this->_inclusions.size()+1);
            _finite_strain_materials.back() = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _finite_strain_materials[i] = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }

            //setup solid material pointer
            _solid_materials.resize(this->_inclusions.size()+1);
            _solid_materials.back() = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_matrix_material);
            for(std::size_t i{0}; i<this->_inclusions.size(); ++i){
                _solid_materials[i] = dynamic_cast<solid_material_base<_T, _Dim, _Container>*>(this->_inclusions[i]->material());
            }

            //setup mixed formulation wrappers
            for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
                const auto& vol_frac{(i == _finite_strain_materials.size() - 1 ? this->_cm : this->_inclusions[i]->volume_fraction())};
                switch (_finite_strain_materials[i]->formulation_type()) {
                case FINITE_STRAIN_FORMULATION::LAGRANGIAN:
                    _reference_formulations.push_back({vol_frac, _solid_materials[i], _finite_strain_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                case FINITE_STRAIN_FORMULATION::MIXED:
                    _mixed_formulations.push_back({vol_frac, _solid_materials[i], this->_strain_concentration_tensors[i]});
                    break;
                }
            }
            solid_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    inline virtual void reinit() override {
        this->_matrix_material->reinit();
        for(auto inclusion : this->_inclusions){
            inclusion->reinit();
        }
        mean_field_composite_solid_base<_T, _Dim, _Container>::reinit();
        update_macro_tangent();
    }

    inline virtual void update_strain()override{
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_strain(): class is not initialized");
        }

        //const auto& F_macro{this->deformation_tensor() + this->_dstrain};
        for(std::size_t i{0}; i<_finite_strain_materials.size(); ++i){
            _finite_strain_materials[i]->deformation_tensor() = tmech::dcontract(this->_strain_concentration_tensors[i], this->deformation_tensor());
        }
    }



    inline virtual void update() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update(): is not initialized");
        }

        update_strain();
        for(auto solid : _solid_materials){
            solid->update();
        }

        update_macro_stress();
        update_macro_tangent();
    }



    inline virtual void update_stress() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_stress(): is not initialized");
        }

        update_strain();

        for(auto solid : _solid_materials){
            solid->update_stress();
        }
        update_macro_stress();
    }



    inline virtual void update_tangent() override {
        if(!this->is_initialized()){
            throw std::runtime_error("small_strain_mean_field_composite::update_tangent(): is not initialized");
        }

        update_strain();
        for(auto solid : _solid_materials){
            solid->update_tangent();
        }
        update_macro_tangent();
    }

    inline virtual solid_material_base<_T, _Dim, _Container>* base_material() override {
        return this;
    }


private:

    constexpr inline auto update_macro_stress(){
        this->_stress.fill(0);
        update_macro_stress_reference();
        update_macro_stress_mixed();
    }

    constexpr inline auto update_macro_stress_reference(){
        for(auto [c, solid, finite, A] : _reference_formulations){
            //const auto& F{finite->deformation_tensor()};
            this->_stress += c*finite->mixed_stress_tensor();//F*solid->stress_tensor();
        }
    }

    constexpr inline auto update_macro_stress_mixed(){
        for(auto [c, solid, A] : _mixed_formulations){
            this->_stress += c*solid->stress_tensor();
        }
    }


    constexpr inline auto update_macro_tangent(){
        this->_C.fill(0);
        update_macro_tangent_reference();
        update_macro_tangent_mixed();
    }

    constexpr inline auto update_macro_tangent_reference(){
        const tmech::eye<value_type, _Dim, 2> I;
        for(auto [c, solid, finite, A] : _reference_formulations){
            //const auto& F{finite->deformation_tensor()};
            this->_C += c*tmech::dcontract(finite->mixed_tangent_tensor(), A);
            //tmech::dcontract(tmech::dcontract(tmech::dcontract(tmech::otimesu(F, I), solid->tangent_tensor()), tmech::otimesu(tmech::trans(F), I))
            //               + tmech::otimesu(I, solid->stress_tensor()), A);
        }
    }

    constexpr inline auto update_macro_tangent_mixed(){
        for(auto [c, solid, A] : _mixed_formulations){
            this->_C += c*tmech::dcontract(solid->tangent_tensor(), A);
        }
    }

    std::vector<finite_strain_solid_material_base<_T, _Dim, _Container>*> _finite_strain_materials;
    std::vector<solid_material_base<_T, _Dim, _Container>*> _solid_materials;
    std::vector<mixed_formulation_wrapper> _mixed_formulations;
    std::vector<reference_formulation_wrapper> _reference_formulations;
};


#endif // FINITE_STRAIN_MEAN_FIELD_COMPOSITE_BONES_H

