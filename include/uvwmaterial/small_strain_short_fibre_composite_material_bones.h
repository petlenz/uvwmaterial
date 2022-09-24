/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef MATERIAL_SMALL_STRAIN_SHORT_FIBRE_COMPOSITE_BONES_H
#define MATERIAL_SMALL_STRAIN_SHORT_FIBRE_COMPOSITE_BONES_H

//short fibre composite
//A_i[p_i] := A[S[p_i]]
//A_m[p_i] = (II - c_i*A_i)/c_m (ausser dilute)
//eps_m = (sum_i w_i A_m_i):eps_macro
//C^{UD}_i := C[v0, 1-v0, C0, C1, p_i]
//Cmacro = sum_i w_i C^{UD}_i:A_i


template <typename _T, std::size_t _Dim, typename _Container>
class small_strain_short_fibre_composite
        : public short_fibre_composite_base<_T, _Dim, _Container>,
          public small_strain_material_base<_T, _Dim, _Container>
{
public:
    using size_type = std::size_t;

    small_strain_short_fibre_composite() {}

    virtual ~small_strain_short_fibre_composite() {}

    virtual inline void init() override {
        if(!small_strain_material_base<_T, _Dim, _Container>::_is_init){
            //init composites
            _small_strain_material.reserve(this->_composites.size());
            for(auto composite : this->_composites){
                _small_strain_material.push_back(dynamic_cast<small_strain_material_base<_T, _Dim, _Container>*>(composite));
                _small_strain_material.back()->init();
            }
            small_strain_material_base<_T, _Dim, _Container>::_is_init = true;
        }
    }

    virtual inline void reinit() override {
        for(auto material : this->_small_strain_material){
            material->reinit();
        }
    }

    virtual inline void update_strain() override {
        for(auto composite : _small_strain_material){
            composite->strain_tensor() = this->_strain;
        }
    }

    virtual inline void update() override {
        update_strain();
        this->_C.fill(0);
        this->_stress.fill(0);
        const auto& weigts{this->weights()};
        for(size_type i{0}; i<weigts.size(); ++i){
            _small_strain_material[i]->update();
            this->_C += weigts[i]*_small_strain_material[i]->tangent_tensor();
            this->_stress += weigts[i]*_small_strain_material[i]->stress_tensor();
        }
    }

    virtual inline void update_tangent() override {
        update_strain();
        this->_C.fill(0);
        const auto& weigts{this->weights()};
        for(size_type i{0}; i<weigts.size(); ++i){
            _small_strain_material[i]->update_tangent();
            this->_C += weigts[i]*_small_strain_material[i]->tangent_tensor();
        }
    }

    virtual inline void update_stress() override {
        update_strain();
        this->_stress.fill(0);
        const auto& weigts{this->weights()};
        for(size_type i{0}; i<weigts.size(); ++i){
            _small_strain_material[i]->update_stress();
            this->_stress += weigts[i]*_small_strain_material[i]->stress_tensor();
        }
    }

    virtual inline solid_material_base<_T, _Dim, _Container>* base_material()override{
        return this;
    }

protected:
    std::vector<small_strain_material_base<_T, _Dim, _Container>*> _small_strain_material;
};






#endif // MATERIAL_SMALL_STRAIN_SHORT_FIBRE_COMPOSITE_BONES_H
