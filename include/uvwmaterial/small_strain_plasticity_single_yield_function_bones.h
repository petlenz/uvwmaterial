/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef SMALL_STRAIN_PLASTICITY_SINGLE_YIELD_FUNCTION_BONES_H
#define SMALL_STRAIN_PLASTICITY_SINGLE_YIELD_FUNCTION_BONES_H

template <typename T, std::size_t Dim, typename Container>
class small_strain_plasticity_single_yield_function :
        public small_strain_material_base<T, Dim, Container>,
        public plastic_material_base<T, Dim>
{
public:
    using value_type = T;

    small_strain_plasticity_single_yield_function():
        small_strain_material_base<T, Dim, Container>(),
        plastic_material_base<T, Dim>(),
        _base_material(nullptr),
        _yield_function(nullptr)
    {
        plastic_material_base<T, Dim>::resize(1);
    }

    small_strain_plasticity_single_yield_function(small_strain_material_base<T, Dim, Container> * __base_material):
        small_strain_material_base<T, Dim, Container>(),
        plastic_material_base<T, Dim>(),
        _base_material(__base_material),
        _yield_function(nullptr)
    {
        plastic_material_base<T, Dim>::resize(1);
    }

    small_strain_plasticity_single_yield_function(small_strain_material_base<T, Dim, Container> * __base_material,
                                                  yield_function_stress_base<value_type, Dim, Container> * __yield_function):
        small_strain_material_base<T, Dim, Container>(),
        plastic_material_base<T, Dim>(),
        _base_material(__base_material),
        _yield_function(__yield_function)
    {
        plastic_material_base<T, Dim>::resize(1);
    }

    constexpr inline void set_base_material(small_strain_material_base<T, Dim, Container> * __base_material){
        _base_material = __base_material;
    }

    constexpr inline void set_yield_function(yield_function_stress_base<value_type, Dim, Container> * __yield_function){
        _yield_function = __yield_function;
    }

    inline virtual void init()override{
        if(!this->_is_init){
            _base_material->init();
            this->_C = _base_material->tangent_tensor();
            this->_is_init = true;
        }
    }

    inline virtual void reinit()override{
        _base_material->reinit();
        this->_C = _base_material->tangent_tensor();
    }

    inline virtual void update()override{
        _base_material->strain_tensor() = this->_strain - this->_inelastic_strain;
        _base_material->update();

        _yield_function->update_equivalent_scalar();

        if(_yield_function->yielding(this->_history[0])){
            const auto G{this->_parameter[0]};
            const auto delta_lambda{_yield_function->solve(this->_history[0])};
            _yield_function->update_stress_dependent_parts();
            const auto& N{_yield_function->normal()};
            this->_inelastic_strain += delta_lambda*N;
            this->_stress = _base_material->stress_tensor() - 2*G*delta_lambda*N;
            this->_C = _base_material->tangent_tensor() - 2*G*tmech::otimes(N, _yield_function->derivative_scalar_wrt_strain()) - 2*G*delta_lambda*_yield_function->derivative_normal_wrt_strain();
            this->_history[0] += delta_lambda;
        }else{
            this->_stress = _base_material->stress_tensor();
            this->_C = _base_material->tangent_tensor();
        }
    }

    inline void update_stress()override{
        _base_material->strain_tensor() = this->_strain - this->_inelastic_strain;
        _base_material->update();

        _yield_function->update_equivalent_scalar();

        if(_yield_function->yielding(this->_history[0])){
            const auto G{this->_parameter[0]};
            const auto delta_lambda{_yield_function->solve(this->_history[0])};
            _yield_function->update_stress_dependent_parts();
            const auto& N{_yield_function->normal()};
            this->_inelastic_strain += delta_lambda*N;
            this->_stress = _base_material->stress_tensor() - 2*G*delta_lambda*N;
            this->_history[0] += delta_lambda;
        }else{
            this->_stress = _base_material->stress_tensor();
        }
    }

    inline void update_tangent()override{
        _base_material->strain_tensor() = this->_strain - this->_inelastic_strain;
        _base_material->update();

        _yield_function->update_equivalent_scalar();

        if(_yield_function->yielding(this->_history[0])){
            const auto G{this->_parameter[0]};
            const auto delta_lambda{_yield_function->solve(this->_history[0])};
            _yield_function->update_stress_dependent_parts();
            const auto& N{_yield_function->normal()};
            this->_inelastic_strain += delta_lambda*N;
            this->_C = _base_material->tangent_tensor() - 2*G*tmech::otimes(N, _yield_function->derivative_scalar_wrt_strain()) - 2*G*delta_lambda*_yield_function->derivative_normal_wrt_strain();
            this->_history[0] += delta_lambda;
        }else{
            this->_C = _base_material->tangent_tensor();
        }
    }

    inline solid_material_base<T, Dim, Container> * base_material()override{
        return _base_material;
    }

protected:
    small_strain_material_base<value_type, Dim, Container>* _base_material;
    yield_function_stress_base<value_type, Dim, Container>* _yield_function;
};



template <typename T, std::size_t Dim, typename Container>
class small_strain_incremental_plasticity_single_yield_function :
        public small_strain_material_base<T, Dim, Container>,
        public plastic_material_base<T, Dim>,
        public incremental_solid_material<T, Dim>
{
public:
    using value_type = T;

    small_strain_incremental_plasticity_single_yield_function():
        small_strain_material_base<T, Dim, Container>(),
        plastic_material_base<T, Dim>(),
        _base_material(nullptr),
        _yield_function(nullptr)
    {
        plastic_material_base<T, Dim>::resize(1);
    }

    small_strain_incremental_plasticity_single_yield_function(small_strain_material_base<T, Dim, Container> * __base_material):
        small_strain_material_base<T, Dim, Container>(),
        plastic_material_base<T, Dim>(),
        _base_material(__base_material),
        _yield_function(nullptr)
    {
        plastic_material_base<T, Dim>::resize(1);
    }

    small_strain_incremental_plasticity_single_yield_function(small_strain_material_base<T, Dim, Container> * __base_material,
                                                              yield_function_stress_base<value_type, Dim, Container> * __yield_function):
        small_strain_material_base<T, Dim, Container>(),
        plastic_material_base<T, Dim>(),
        _base_material(__base_material),
        _yield_function(__yield_function)
    {
        plastic_material_base<T, Dim>::resize(1);
    }

    constexpr inline void set_base_material(small_strain_material_base<T, Dim, Container> * __base_material){
        _base_material = __base_material;
    }

    constexpr inline void set_yield_function(yield_function_stress_base<value_type, Dim, Container> * __yield_function){
        _yield_function = __yield_function;
    }

    inline void init()override{
        if(!this->_is_init){
            _incremental_base_material = dynamic_cast<incremental_solid_material<value_type, Dim>*>(_base_material);
            if(!_incremental_base_material){
                throw std::runtime_error("small_strain_incremental_plasticity_single_yield_function::init(): base material is not an incremental materail");
            }
            _base_material->init();
            this->_C = _base_material->tangent_tensor();
            this->_is_init = true;
        }
    }

    inline void reinit()override{
        _base_material->reinit();
        this->_C = _base_material->tangent_tensor();
    }

    inline void update()override{
        _incremental_base_material->dstrain_tensor() = this->dstrain_tensor();
        _base_material->stress_tensor() = this->stress_tensor();
        _base_material->update();

        _yield_function->update_equivalent_scalar();

        if(_yield_function->yielding(this->_history[0])){
            const auto G{this->_parameter[0]};
            const auto delta_lambda{_yield_function->solve(this->_history[0])};
            _yield_function->update_stress_dependent_parts();
            const auto& N{_yield_function->normal()};
            this->_inelastic_strain += delta_lambda*N;
            this->_stress = _base_material->stress_tensor() - 2*G*delta_lambda*N;
            this->_C = _base_material->tangent_tensor() - 2*G*tmech::otimes(N, _yield_function->derivative_scalar_wrt_strain()) - 2*G*delta_lambda*_yield_function->derivative_normal_wrt_strain();
            this->_history[0] += delta_lambda;
        }else{
            this->_stress = _base_material->stress_tensor();
            this->_C = _base_material->tangent_tensor();
        }
    }

    inline void update_stress()override{
        _incremental_base_material->dstrain_tensor() = this->dstrain_tensor();
        _base_material->stress_tensor() = this->stress_tensor();
        _base_material->update();

        _yield_function->update_equivalent_scalar();

        if(_yield_function->yielding(this->_history[0])){
            const auto G{this->_parameter[0]};
            const auto delta_lambda{_yield_function->solve(this->_history[0])};
            _yield_function->update_stress_dependent_parts();
            const auto& N{_yield_function->normal()};
            this->_inelastic_strain += delta_lambda*N;
            this->_stress = _base_material->stress_tensor() - 2*G*delta_lambda*N;
            this->_history[0] += delta_lambda;
        }else{
            this->_stress = _base_material->stress_tensor();
        }
    }

    inline void update_tangent()override{
        _incremental_base_material->dstrain_tensor() = this->dstrain_tensor();
        _base_material->stress_tensor() = this->stress_tensor();
        _base_material->update();

        _yield_function->update_equivalent_scalar();

        if(_yield_function->yielding(this->_history[0])){
            const auto G{this->_parameter[0]};
            const auto delta_lambda{_yield_function->solve(this->_history[0])};
            _yield_function->update_stress_dependent_parts();
            const auto& N{_yield_function->normal()};
            this->_inelastic_strain += delta_lambda*N;
            this->_C = _base_material->tangent_tensor() - 2*G*tmech::otimes(N, _yield_function->derivative_scalar_wrt_strain()) - 2*G*delta_lambda*_yield_function->derivative_normal_wrt_strain();
            this->_history[0] += delta_lambda;
        }else{
            this->_C = _base_material->tangent_tensor();
        }
    }

    inline solid_material_base<T, Dim, Container> * base_material()override{
        return _base_material;
    }

protected:
    small_strain_material_base<value_type, Dim, Container>* _base_material;
    yield_function_stress_base<value_type, Dim, Container>* _yield_function;
    incremental_solid_material<value_type, Dim>* _incremental_base_material;;
};
#endif // SMALL_STRAIN_PLASTICITY_SINGLE_YIELD_FUNCTION_BONES_H
