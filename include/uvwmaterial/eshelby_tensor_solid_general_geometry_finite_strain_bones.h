/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
//#ifndef ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_FINITE_STRAIN_BONES_H
//#define ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_FINITE_STRAIN_BONES_H

//template <typename _T, std::size_t _Dim, typename _Container>
//class eshelby_tensor_solid_general_geometry_finite_strain :
//        public eshelby_tensor_solid_base<_T, _Dim>,
//        //public eshelby_tensor_solid_general_geometry<_T, _Dim>,
//        public eshelby_tensor_finite_strain_base<_T, _Dim, _Container>
//{
//public:
//    using value_type = _T;
//    using eshelby_nonlinear_base = eshelby_tensor_finite_strain_base<_T, _Dim, _Container>;
//    using eshelby_base           = eshelby_tensor_solid_base<_T, _Dim>;

//    eshelby_tensor_solid_general_geometry_finite_strain()
//    {}

//    virtual ~eshelby_tensor_solid_general_geometry_finite_strain()
//    {}

//    virtual void init() override{
//        if(!this->_is_init){
//            //eshelby_base::init();
//            eshelby_nonlinear_base::init();
//            eshelby_nonlinear_base::_parameter.resize(_Dim+1);
//            reinit();
////            this->update_isotropic_matrix();
////            eshelby_base::_parameter[0] = eshelby_nonlinear_base::_parameter[3];
////            eshelby_nonlinear_base::_parameter[0] = eshelby_base::_parameter[1];
////            eshelby_nonlinear_base::_parameter[1] = eshelby_base::_parameter[2];
////            eshelby_nonlinear_base::_parameter[2] = eshelby_base::_parameter[3];
////            eshelby_nonlinear_base::init();
//            //_S0 = this->_S;
//            this->_is_init = true;
//        }
//    }

//    virtual void reinit() override {
//        this->update_geometry();
//        //this->update_isotropic_matrix();
//        //eshelby_base::_parameter[0] = eshelby_nonlinear_base::_parameter[3];
//        //eshelby_base::_parameter[1] = eshelby_nonlinear_base::_parameter[0];
//        //eshelby_base::_parameter[2] = eshelby_nonlinear_base::_parameter[1];
//        //eshelby_base::_parameter[3] = eshelby_nonlinear_base::_parameter[2];
//        //eshelby_tensor_solid_general_geometry<_T, _Dim>::reinit();

//        this->num_evaluate();
//        eshelby_base::_S = eshelby_nonlinear_base::_S;

//        //this->_S = (tmech::dcontract(tmech::otimesu(this->_R,this->_R),tmech::dcontract(this->_S,tmech::otimesu(tmech::trans(this->_R),tmech::trans(this->_R)))));
////        auto matrix_material = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(this->_base_material);
////        //this->_C = matrix_material->spatial_tangent_tensor();
////        tmech::tensor<_T,_Dim,4> A = matrix_material->mixed_tangent_tensor(), Ainv;
//////        Ainv = tmech::invf(A);
//////        if(std::isnan(Ainv(0,0,0,0))){
//////            Ainv = tmech::inv(A);
//////        }
////        auto is_minor_sym = [](auto const& A){
////            constexpr auto Dim{A.dimension()};
////            for(std::size_t i{0}; i<Dim; ++i){
////                for(std::size_t j{0}; j<Dim; ++j){
////                    for(std::size_t k{0}; k<Dim; ++k){
////                        for(std::size_t l{0}; l<Dim; ++l){
////                            if(std::abs(A(i,j,k,l) - A(j,i,k,l)) > std::numeric_limits<value_type>::epsilon()*std::max((_T)1.0, std::max(A(i,j,k,l), A(j,i,k,l)))*10.0 ||
////                                    std::abs(A(i,j,k,l) - A(i,j,l,k)) > std::numeric_limits<value_type>::epsilon()*std::max((_T)1.0, std::max(A(i,j,k,l), A(i,j,l,k)))*10.0 ){
////                                //std::cout<<std::abs(A(i,j,k,l) - A(j,i,k,l))<<" "<<std::abs(A(i,j,k,l) - A(i,j,l,k))<<" "<<std::numeric_limits<value_type>::epsilon()*std::max(A(i,j,k,l), A(j,i,k,l))<<std::endl;
////                                return false;
////                            }
////                        }
////                    }
////                }
////            }
////            return true;
////        };
////        if(is_minor_sym(A)){
////            Ainv = tmech::inv(A);
////        }else{
////            Ainv = tmech::invf(A);
////        }
////        this->_S = tmech::eval(tmech::dcontract(tmech::dcontract(Ainv,tmech::dcontract(this->_C,tmech::dcontract(this->_S,tmech::inv(this->_C)))),A));
////        std::cout<<this->_S<<std::endl;
//    }

////    virtual void update() override {
//////        this->update_geometry();
//////        this->update_isotropic_matrix();
//////        eshelby_base::_parameter[0] = eshelby_nonlinear_base::_parameter[0];
//////        eshelby_base::_is_init = false;
//////        eshelby_base::init();
////    }

//private:
//    //tmech::tensor<_T,_Dim,4> _S0;
//};



//#endif // ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_FINITE_STRAIN_BONES_H
