/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_FINITE_STRAIN_BASE_BONES_H
#define ESHELBY_TENSOR_FINITE_STRAIN_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class eshelby_tensor_finite_strain_base :
        public eshelby_tensor_solid_base<_T, _Dim>
{
public:
    using value_type = _T;
    using integration_data = SPHERE_LEBEDEV_RULE<value_type>;

    eshelby_tensor_finite_strain_base():
        _A0(),
        _R(),//tmech::eye<_T,_Dim,2>()),
        _base_material(nullptr)
    {}

    virtual ~eshelby_tensor_finite_strain_base(){}

    //virtual void update() = 0;

    constexpr inline void init(){
        if(!_base_material){
            throw std::runtime_error("eshelby_tensor_finite_strain_base::init() base material is not set");
        }
        _parameter_update.resize(_Dim);
        update_geometry();
        num_evaluate();
    }

    constexpr inline void reinit(){
        init();
    }

    constexpr inline auto update_geometry(){
        //use deformation gradient from _material_base
        //std::cout<<"start"<<std::endl;
        for(std::size_t i{0}; i<_Dim; ++i){
            _A0(i,i) = 1.0/(this->_parameter[i]*this->_parameter[i]);
            //std::cout<<this->_parameter[i]<<"\n";
        }
        //std::cout<<std::endl;
        const auto& F = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(_base_material)->deformation_tensor();
        const auto Finv{tmech::inv(F)};
        auto Aeig{tmech::eigen_decomposition(tmech::trans(Finv)*_A0*Finv)};
        Aeig.decompose();
        //std::cout<<"Start geometry"<<std::endl;
        //std::cout<<this->_parameter[0]<<" "<<this->_parameter[1]<<" "<<this->_parameter[2]<<std::endl;
        this->_parameter_update[0] = std::sqrt(1/Aeig.eigenvalues()[0]);//[_index[0]]);
        this->_parameter_update[1] = std::sqrt(1/Aeig.eigenvalues()[1]);//[_index[1]]);
        this->_parameter_update[2] = std::sqrt(1/Aeig.eigenvalues()[2]);//[_index[2]]);
        //        if(this->_parameter[0] != 3){
        //            std::cout<<"H"<<std::endl;
        //        }
        //std::cout<<Aeig.permutation()[0]<<" "<<Aeig.permutation()[1]<<" "<<Aeig.permutation()[2]<<std::endl;
        //std::cout<<std::sqrt(1/Aeig.eigenvalues()[0])<<" "<<std::sqrt(1/Aeig.eigenvalues()[1])<<" "<<std::sqrt(1/Aeig.eigenvalues()[2])<<std::endl;
        //std::cout<<this->_parameter_update[0]<<" "<<this->_parameter_update[1]<<" "<<this->_parameter_update[2]<<std::endl;





        for(std::size_t i{0}; i<_Dim; ++i){
            for(std::size_t j{0}; j<_Dim; ++j){
                _R(i,j) = Aeig.eigenvectors()[j](i);
            }
        }

        if constexpr(_Dim == 3){
            tmech::adaptor<_T,_Dim,1,tmech::full<_Dim>> geo(&this->_parameter_update[0]);
            _geometry = geo;
        }else{
            throw;
        }


        //tmech::tensor<value_type,_Dim,1> eigval{Aeig.eigenvalues()[0],Aeig.eigenvalues()[1],Aeig.eigenvalues()[2]};
        //std::cout<<_R<<std::endl;
        //eigval = tmech::eval(_R*eigval);//*tmech::trans(_R));
        //std::cout<<std::sqrt(1/eigval(0))<<" "<<std::sqrt(1/eigval(1))<<" "<<std::sqrt(1/eigval(2))<<std::endl;
        //std::cout<<std::sqrt(1/std::abs(eigval(0)))<<" "<<std::sqrt(1/std::abs(eigval(1)))<<" "<<std::sqrt(1/std::abs(eigval(2)))<<std::endl;

        if(_R(0,2) != 0 && 0 != _R(2,0)){
            //std::cout<<std::atan(_R_temp(1,2)/_R_temp(0,2))<<" "<<std::acos(_R_temp(2,2))<<" "<<std::atan(-_R_temp(2,1)/_R_temp(2,0))<<std::endl;
        }
        //        if(tmech::norm(_R) != 0){
        //            _R = tmech::eval(_R_temp*_R);
        //        }else{
        //            _R = _R_temp;
        //        }
        //_R += tmech::eval(_R_temp*_R);
    }

    constexpr inline auto update_isotropic_matrix(){
        const tmech::eye<value_type,_Dim,2> I;
        //fourth order symmetric idenity tensor
        const auto IIsym = 0.5*(tmech::otimesu(I,I)+tmech::otimesl(I,I));
        //fourth order identity tensor
        const auto II = tmech::otimes(I,I);
        const auto IIvol = II/_Dim;
        const auto IIdev = IIsym-IIvol;
        auto matrix_material = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(_base_material);
        _C = matrix_material->spatial_tangent_tensor();
        const auto K{tmech::ddcontract(IIvol, _C)/3};
        const auto G{tmech::ddcontract(IIdev, _C)/10};
        std::cout<<(_C(3,3,3,3) + _C(4,4,4,4) + _C(5,5,5,5))/3.<<std::endl;
        this->_parameter[_Dim] = (3*K-2*G)/(2*(3*K+G));
        std::cout<<"K = "<<K<<" G = "<<G<<" "<<this->_parameter[_Dim]<<std::endl;
        if(this->_parameter[_Dim] < 0 || this->_parameter[_Dim] > 0.5){
            //std::cout<<"K = "<<_K<<" G = "<<G<<std::endl;
            throw std::runtime_error("eshelby_tensor_finite_strain_base::update_isotropic_matrix() negative parameter nu = " + std::to_string(this->_parameter[_Dim]));
        }
    }


    constexpr inline auto num_evaluate(){
        auto matrix_material = dynamic_cast<finite_strain_solid_material_base<_T, _Dim, _Container>*>(_base_material);
        tmech::tensor<_T,_Dim,4> F = (tmech::dcontract(tmech::otimesu(this->_R,this->_R), tmech::dcontract(matrix_material->mixed_tangent_tensor(), tmech::otimesu(tmech::trans(this->_R),tmech::trans(this->_R)))));
        _P.fill(0);
        const auto size{integration_data::get_size(11)};
        const value_type* phi{integration_data::get_phi(11)};
        const value_type* theta{integration_data::get_theta(11)};
        const value_type* weights{integration_data::get_weight(11)};
        const auto a1{this->_parameter_update[0]}, a2{this->_parameter_update[1]}, a3{this->_parameter_update[2]};
        const tmech::tensor<_T,_Dim,4> Ft{tmech::basis_change<tmech::sequence<1,3,2,4>>(F)};
        for(std::size_t ii{0}; ii<size; ++ii){
            const auto _theta{theta[ii]/180*M_PI};
            const auto _phi{phi[ii]/180*M_PI};
            const tmech::tensor<_T,_Dim,1> x{std::cos(_theta)*std::sin(_phi)/a1, std::sin(_theta)*std::sin(_phi)/a2, std::cos(_phi)/a3};
            const tmech::tensor<_T,_Dim,2> Ainv{tmech::inv(tmech::dcontract(tmech::otimes(x,x), Ft))};
            const auto w = weights[ii];
            for(std::size_t i{0}; i<_Dim; ++i){
                for(std::size_t j{0}; j<_Dim; ++j){
                    for(std::size_t k{0}; k<_Dim; ++k){
                        for(std::size_t l{0}; l<_Dim; ++l){
                            _P(i,j,k,l) += Ainv(i,k)*x(j)*x(l)*w;
                        }
                    }
                }
            }
        }


        //        std::size_t Ntheta{64};
        //        const _T dphi{2*M_PI/(2*Ntheta)};
        //        const _T dtheta{M_PI/(Ntheta)};
        //        for(std::size_t itheta{0}; itheta<Ntheta; ++itheta){
        //            _T theta{itheta*M_PI/Ntheta};
        //            //std::cout<<"theta "<<theta<<std::endl;
        //            for(std::size_t iphi{0}; iphi<Ntheta*2; ++iphi){
        //                _T phi{2*iphi*M_PI/(2*Ntheta)};
        //                //std::cout<<"phi "<<phi<<std::endl;
        //                tmech::tensor<_T,_Dim,1> x{std::sin(theta)*std::cos(phi)/a1, std::sin(theta)*std::sin(phi)/a2, std::cos(theta)/a3};
        //                const tmech::tensor<_T,_Dim,2> invA{tmech::inv(x*tmech::transl(F)*x)};
        //                const auto w = 4*M_PI*dtheta*dphi;
        //                integrate(F,invA,x,w);
        //            }
        //        }
        this->_S = (tmech::dcontract(tmech::otimesu(tmech::trans(this->_R), tmech::trans(this->_R)),tmech::dcontract(tmech::dcontract(_P,F),tmech::otimesu((this->_R),(this->_R)))));
    }

    //    template<typename _Ainv, typename _x>
    //    constexpr inline auto integrate(_Ainv const& Ainv, _x const& x, _T w){
    //        for(std::size_t i{0}; i<_Dim; ++i){
    //            for(std::size_t j{0}; j<_Dim; ++j){
    //                for(std::size_t k{0}; k<_Dim; ++k){
    //                    for(std::size_t l{0}; l<_Dim; ++l){
    //                        _P(i,j,k,l) += (Ainv(j,k)*x(i)*x(l) + Ainv(i,k)*x(j)*x(l) + Ainv(j,l)*x(i)*x(k) + Ainv(i,l)*x(j)*x(k))*w;
    //                    }
    //                }
    //            }
    //        }
    //    }

    constexpr inline auto set_base_material(material_base<_T, _Dim, _Container>* __base_material){
        _base_material = __base_material;
    }

    constexpr inline auto const& updated_parameter()const{
        return _parameter_update;
    }

    constexpr inline auto& updated_parameter(){
        return _parameter_update;
    }

    constexpr inline auto const& geometry()const{
        return _geometry;
    }

    constexpr inline auto const& rotation_tensor()const{
        return _R;
    }

    constexpr inline auto& rotation_tensor(){
        return _R;
    }

protected:
    _T _K0;
    tmech::tensor<_T, _Dim, 2> _A0;
    tmech::tensor<_T, _Dim, 2> _R;
    tmech::tensor<_T, _Dim, 4> _C;
    tmech::tensor<_T, _Dim, 4> _S;
    tmech::tensor<_T, _Dim, 4> _P;
    tmech::tensor<_T, _Dim, 1> _geometry;
    material_base<_T, _Dim, _Container>* _base_material;
    _Container _parameter_update;
};

#endif // ESHELBY_TENSOR_FINITE_STRAIN_BASE_BONES_H
