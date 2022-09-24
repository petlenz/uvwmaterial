/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef COMPOSITE_MATERIAL_EIGENSTRAIN_BASE_BONES_H
#define COMPOSITE_MATERIAL_EIGENSTRAIN_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class composite_eigenstrain_material_base :
        public eigenstrain_material_base<_T, _Dim>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    composite_eigenstrain_material_base(material_base<_T, _Dim, _Container> * __material);

    composite_eigenstrain_material_base(material_base<_T, _Dim, _Container> * __material, size_type const __number_of_inclusions, size_type const __number_of_eigenstrains);

    virtual ~composite_eigenstrain_material_base(){}

    constexpr inline auto init(){
        const auto material_composite{dynamic_cast<composite_material_base<_T, _Dim, _Container>*>(_material)};
        const auto& C{_material->tangent_tensor()};
        const auto& strain_concentration_tensors{material_composite->get_strain_concentration_tensors()};
        const auto& A_m{strain_concentration_tensors.back()};
        const auto c_m{material_composite->get_matrix_volume_fraction()};
        const auto& C_m{material_composite->get_matrix_material()->tangent_tensor()};
        const auto I{tmech::eye<value_type, _Dim, 2>()};
        const auto IIsym{(tmech::otimesu(I,I) + (tmech::otimesl(I,I)))*0.5};
        const auto& inclusions{material_composite->get_inclusions()};
        this->eigenstrain_concentration_tensors.resize(inclusions.size()+1);
        for(size_type i{0}; i<inclusions.size(); ++i){
            this->eigenstrain_concentration_tensors[i].resize(inclusions.size()+1);
            const auto& A_i{strain_concentration_tensors[i]};
            const auto c_i{inclusions[i]->get_volume_fraction()};
            const auto& C_i{inclusions[i]->get_material()->get_tangent_tensor()};
            //D_{ir} = (IIsym - A_i):inv(C_i - C):(delta_ir*IIsym - c_i*trans(A_i)):C_i
            const tmech::tensor<value_type, _Dim, 4> D_ir{tmech::dcontract((IIsym - A_i), tmech::dcontract(tmech::inv(C_i - C), tmech::dcontract( - c_i*tmech::trans(A_i), C_i)))};
            for(size_type r{0}; r<inclusions.size(); ++r){
                if(i == r){
                    this->eigenstrain_concentration_tensors[i][r] = tmech::dcontract((IIsym - A_i), tmech::dcontract(tmech::inv(C_i - C), tmech::dcontract((IIsym - c_i*tmech::trans(A_i)), C_i)));
                }else{
                    this->eigenstrain_concentration_tensors[i][r] = D_ir;
                }
            }
            //matrix part
            this->eigenstrain_concentration_tensors[i].back() = D_ir;
        }
        //matrix phase

        this->eigenstrain_concentration_tensors.back().resize(inclusions.size()+1);
        const tmech::tensor<value_type, _Dim, 4> D_ir{tmech::dcontract((IIsym - A_m), tmech::dcontract(tmech::inv(C_m - C), tmech::dcontract( - c_m*tmech::trans(A_m), C_m)))};
        for(size_type r{0}; r<inclusions.size(); ++r){
            this->eigenstrain_concentration_tensors.back()[r] = D_ir;
        }
        this->eigenstrain_concentration_tensors.back().back() = tmech::dcontract((IIsym - A_m), tmech::dcontract(inv(C_m - C), tmech::dcontract((IIsym - c_m*tmech::trans(A_m)), C_m)));

        //stress concentration tensors
        const tmech::tensor<value_type, _Dim, 4> Cinv{tmech::inv(C)};
        this->stress_concentration_tensors.resize(inclusions.size()+1);
        for(size_type i{0}; i<inclusions.size(); ++i){
            const auto& C_i{inclusions[i]->get_material()->get_tangent_tensor()};
            this->stress_concentration_tensors[i] = tmech::dcontract(C_i, tmech::dcontract(strain_concentration_tensors[i], Cinv));
        }
        this->stress_concentration_tensors.back() = tmech::dcontract(C_m, tmech::dcontract(strain_concentration_tensors.back(), Cinv));
    }

//    inline void update_stress(){
//        update_strain();
//        const auto material_composite{dynamic_cast<composite_material_base<T, Dim, Container>*>(material)};
//        auto small_strain_composite_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material)};
//        const auto& inclusions{material_composite->get_inclusions()};
//        auto matrix_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material_composite->get_matrix_material())};

//        matrix_material->update();
//        const_cast<tmech::tensor<value_type, Dim, 2>&>(small_strain_composite_material->get_stress_tensor()) = matrix_material->get_stress_tensor()*material_composite->get_matrix_volume_fraction();
//        for(size_type i{0}; i<inclusions.size(); ++i){
//            auto small_strain_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(inclusions[i]->get_material())};
//            small_strain_material->update();
//            const_cast<tmech::tensor<value_type, Dim, 2>&>(small_strain_composite_material->get_stress_tensor()) += small_strain_material->get_stress_tensor()*inclusions[i]->get_volume_fraction();
//        }
//    }

//    inline void update_tangent(){
//        update_strain();
//        const auto material_composite{dynamic_cast<composite_material_base<T, Dim, Container>*>(material)};
//        auto small_strain_composite_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material)};
//        const auto& inclusions{material_composite->get_inclusions()};
//        auto matrix_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material_composite->get_matrix_material())};
//        const auto& strain_concentration_tensors{material_composite->get_strain_concentration_tensors()};

//        matrix_material->update();
//        const_cast<tmech::tensor<value_type, Dim, 4>&>(small_strain_composite_material->get_tangent_tensor()) = material_composite->get_matrix_volume_fraction()*dcontract(matrix_material->get_tangent_tensor(), strain_concentration_tensors.back());
//        for(size_type i{0}; i<inclusions.size(); ++i){
//            auto small_strain_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(inclusions[i]->get_material())};
//            small_strain_material->update();
//            const_cast<tmech::tensor<value_type, Dim, 4>&>(small_strain_composite_material->get_tangent_tensor()) += inclusions[i]->get_volume_fraction()*dcontract(small_strain_material->get_tangent_tensor(), strain_concentration_tensors[i]);
//        }
//    }

//    inline void update(){
//        update_strain();
//        const auto material_composite{dynamic_cast<composite_material_base<T, Dim, Container>*>(material)};
//        auto small_strain_composite_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material)};
//        const auto& inclusions{material_composite->get_inclusions()};
//        auto matrix_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material_composite->get_matrix_material())};
//        const auto& strain_concentration_tensors{material_composite->get_strain_concentration_tensors()};


//        matrix_material->update();
//        const_cast<tmech::tensor<value_type, Dim, 2>&>(small_strain_composite_material->get_stress_tensor()) = matrix_material->get_stress_tensor()*material_composite->get_matrix_volume_fraction();
//        const_cast<tmech::tensor<value_type, Dim, 4>&>(small_strain_composite_material->get_tangent_tensor())  = material_composite->get_matrix_volume_fraction()*dcontract(matrix_material->get_tangent_tensor(), strain_concentration_tensors.back());

//        for(size_type i{0}; i<inclusions.size(); ++i){
//            auto small_strain_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(inclusions[i]->get_material())};
//            small_strain_material->update();
//            const_cast<tmech::tensor<value_type, Dim, 2>&>(small_strain_composite_material->get_stress_tensor()) += small_strain_material->get_stress_tensor()*inclusions[i]->get_volume_fraction();
//            const_cast<tmech::tensor<value_type, Dim, 4>&>(small_strain_composite_material->get_tangent_tensor())  += inclusions[i]->get_volume_fraction()*dcontract(small_strain_material->get_tangent_tensor(), strain_concentration_tensors[i]);
//        }

//        update_eigenstrain();
//    }


//    inline void update_strain(){
//        const auto material_composite{dynamic_cast<composite_material_base<T, Dim, Container>*>(material)};
//        auto small_strain_composite_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material)};
//        const auto& inclusions{material_composite->get_inclusions()};

//        const auto& eps_macro{small_strain_composite_material->get_strain_tensor()};
//        auto matrix_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material_composite->get_matrix_material())};
//        matrix_material->set_strain_tensor(dcontract(material_composite->get_strain_concentration_tensors().back(), eps_macro) + determine_eigenstrain_tensor(inclusions.size()));
//        for(size_type i{0}; i<inclusions.size(); ++i){
//            auto small_strain_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(inclusions[i]->get_material())};
//            small_strain_material->set_strain_tensor(dcontract(material_composite->get_strain_concentration_tensors()[i], eps_macro)  + determine_eigenstrain_tensor(i));
//        }
//    }




//private:
//    inline void update_eigenstrain(){
//        const auto material_composite{dynamic_cast<composite_material_base<T, Dim, Container>*>(material)};
//        auto small_strain_composite_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material)};
//        auto matrix_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material_composite->get_matrix_material())};
//        const auto& strain_concentration_tensors{material_composite->get_strain_concentration_tensors()};
//        const auto& inclusions{material_composite->get_inclusions()};

//        auto eigenstrain_material{dynamic_cast<uvwmat::eigenstrain_material_base<value_type, Dim>*>(matrix_material->get_base_material())};

//        if(eigenstrain_material){
//            this->eigenstrain_tensor = material_composite->get_matrix_volume_fraction()*dcontract(tmech::trans(this->stress_concentration_tensors.back()), eigenstrain_material->get_eigenstrain_tensor());
////            for(size_type i{0}; i<eigenstrain_material->get_eigenstrain_tensors().size(); ++i){
////            }
//        }

//        for(size_type i{0}; i<inclusions.size(); ++i){
//            auto eigenstrain_material{dynamic_cast<uvwmat::eigenstrain_material_base<value_type, Dim>*>(inclusions[i]->get_material()->get_base_material())};
//            if(eigenstrain_material){
//                this->eigenstrain_tensor += inclusions[i]->get_volume_fraction()*dcontract(tmech::trans(this->stress_concentration_tensors[i]), eigenstrain_material->get_eigenstrain_tensor());
////                for(size_type j{0}; j<eigenstrain_material->get_eigenstrain_tensors().size(); ++j){
////                    this->eigenstrain_tensors[j] += inclusions[i]->get_volume_fraction()*dcontract(transD(this->stress_concentration_tensors[i]), eigenstrain_material->get_eigenstrain_tensors()[j]);
////                }
//            }
//        }
//    }

//    constexpr inline auto determine_eigenstrain_tensor(size_type const i)const {
//        const auto material_composite{dynamic_cast<composite_material_base<T, Dim, Container>*>(material)};
//        auto matrix_material{dynamic_cast<small_strain_material_base<T, Dim, Container>*>(material_composite->get_matrix_material())};
//        tmech::tensor<value_type, Dim, 2> eig_strain;

//        auto eigenstrain_material{dynamic_cast<uvwmat::eigenstrain_material_base<value_type, Dim>*>(matrix_material->get_base_material())};
//        if(eigenstrain_material){
//            eig_strain += dcontract(this->eigenstrain_concentration_tensors[i].back(), eigenstrain_material->get_eigenstrain_tensor());
//        }

//        const auto& inclusions{material_composite->get_inclusions()};
//        for(size_type j{0}; j<inclusions.size(); ++j){
//            //check if material is eigenstrain
//            auto eigenstrain_material{dynamic_cast<uvwmat::eigenstrain_material_base<value_type, Dim>*>(inclusions[j]->get_material()->get_base_material())};
//            if(eigenstrain_material){
//                eig_strain += dcontract(this->eigenstrain_concentration_tensors[i][j], eigenstrain_material->get_eigenstrain_tensor());
////                for(const auto& tensor : eigenstrain_material->get_eigenstrain_tensors()){
////                    eig_strain += dcontract(this->eigenstrain_concentration_tensors[i][j], tensor);
////                }
//            }
//        }

//        return eig_strain;
//    }

private:
    material_base<_T, _Dim, _Container> * _material;

protected:
    std::vector<tmech::tensor<value_type, _Dim, 4>> _stress_concentration_tensors;
    std::vector<std::vector<tmech::tensor<value_type, _Dim, 4>>> _eigenstrain_concentration_tensors;
};

#endif // COMPOSITE_MATERIAL_EIGENSTRAIN_BASE_BONES_H
