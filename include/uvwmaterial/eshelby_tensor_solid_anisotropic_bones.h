/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_ANISOTROPIC_BONES_H
#define ESHELBY_TENSOR_SOLID_ANISOTROPIC_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class eshelby_tensor_solid_anisotropic :
        public eshelby_tensor<_T, _Dim>
{
public:
    eshelby_tensor_solid_anisotropic():
        _material_base(nullptr)
    {}

    inline void init(){
        throw std::runtime_error("eshelby_tensor_anisotropic not implemented");
    }

    inline virtual void reinit(){

    }
private:
    material_base<_T, _Dim, _Container>* _material_base;
};

//namespace details {

//struct eshelby_tensor_fiber_roots_weights
//{
//    //number of inner integration points
//    static constexpr std::size_t number_inner{16};
//    //number of outer integration points
//    static constexpr std::size_t number_outer{2};

//    static constexpr double roots_inner[]{-0.9894009349916495,
//                                   -0.9445750230732327,
//                                   -0.8656312023878319,
//                                   -0.7554044083550033,
//                                   -0.6178762444026443,
//                                   -0.4580167776572275,
//                                   -0.2816035507792588,
//                                   -0.0950125098376372,
//                                   0.0950125098376372,
//                                   0.2816035507792590,
//                                   0.4580167776572271,
//                                   0.6178762444026440,
//                                   0.7554044083550031,
//                                   0.8656312023878318,
//                                   0.9445750230732326,
//                                   0.9894009349916498};

//    static constexpr double  weights_inner[]{0.0271524594117546,
//                                     0.0622535239386470,
//                                     0.0951585116824926,
//                                     0.1246289712555342,
//                                     0.1495959888165764,
//                                     0.1691565193950028,
//                                     0.1826034150449246,
//                                     0.1894506104550683,
//                                     0.1894506104550689,
//                                     0.1826034150449239,
//                                     0.1691565193950028,
//                                     0.1495959888165774,
//                                     0.1246289712555348,
//                                     0.0951585116824925,
//                                     0.0622535239386481,
//                                     0.0271524594117544};

//    static constexpr double roots_outer[]{-0.57735,
//                                           0.57735};

//    static constexpr double  weights_outer[]{1.0,
//                                             1.0};
//};

//}

//template<typename T>
//class eshelby_tensor
//{
//public:
//    using value_type = T;
//    using size_type = std::size_t;

//    eshelby_tensor() {}

//    template<typename Tensor, typename Geometriy>
//    constexpr inline auto set_up(Tensor const& elasticity_tensor, Geometriy const& a){
//        data_.fill(0);
//        constexpr auto roots_weights{details::eshelby_tensor_fiber_roots_weights()};

//        eshelby_index_notation_fast(roots_weights,a,elasticity_tensor);
//    }

//    template<typename Geometry, typename IntegrationPoints, typename Tensor>
//    constexpr inline auto eshelby_index_notation_fast(IntegrationPoints const& points, Geometry const& a, Tensor const& L){
//        std::array<value_type,36> sum{0};
//        uvw::math::matrix<value_type,3,3> K;
//        constexpr value_type eigth{1.0/8.0};

//        for(size_type M{0};M<points.number_outer;++M){
//            const auto Zeta3 = points.roots_outer[M];
//            const auto W_M = points.weights_outer[M];
//            for(size_type N{0};N<points.number_inner;++N){
//                const auto omega_q = M_PI*(points.roots_inner[N]+1);
//                const auto W_N = points.weights_inner[N];
//                const std::array<value_type,3> Xi{std::sqrt(1.0-Zeta3*Zeta3)*std::cos(omega_q)/a[0],std::sqrt(1.0-Zeta3*Zeta3)*std::sin(omega_q)/a[1],Zeta3/a[2]};//Xi(a,omega_q,Zeta3)};
//                compute_K(K,L,Xi);
//                value_type const D{compute_D(K)};
//                const value_type weights{W_M*W_N};
//                for(size_type m{0};m<3;++m){
//                    for(size_type n{0};n<3;++n){
//                        //first row
//                        {
//                            constexpr size_type i{0},j{0},k{0},l{0};
//                            sum[0] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{0},k{1},l{1};
//                            sum[1] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{0},k{2},l{2};
//                            sum[2] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{0},k{1},l{2};
//                            sum[3] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{0},k{0},l{2};
//                            sum[4] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{0},k{0},l{1};
//                            sum[5] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }

//                        //second row
//                        {
//                            constexpr size_type i{1},j{1},k{0},l{0};
//                            sum[6] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{1},k{1},l{1};
//                            sum[7] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{1},k{2},l{2};
//                            sum[8] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{1},k{1},l{2};
//                            sum[9] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{1},k{0},l{2};
//                            sum[10] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{1},k{0},l{1};
//                            sum[11] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }

//                        //third rows
//                        {
//                            constexpr size_type i{2},j{2},k{0},l{0};
//                            sum[12] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{2},j{2},k{1},l{1};
//                            sum[13] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{2},j{2},k{2},l{2};
//                            sum[14] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{2},j{2},k{1},l{2};
//                            sum[15] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{2},j{2},k{0},l{2};
//                            sum[16] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{2},j{2},k{0},l{1};
//                            sum[17] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }

//                        //fourth row
//                        {
//                            constexpr size_type i{1},j{2},k{0},l{0};
//                            sum[18] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{2},k{1},l{1};
//                            sum[19] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{2},k{2},l{2};
//                            sum[20] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{2},k{1},l{2};
//                            sum[21] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{2},k{0},l{2};
//                            sum[22] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{1},j{2},k{0},l{1};
//                            sum[23] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }

//                        //fifth row
//                        {
//                            constexpr size_type i{0},j{2},k{0},l{0};
//                            sum[24] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{2},k{1},l{1};
//                            sum[25] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{2},k{2},l{2};
//                            sum[26] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{2},k{1},l{2};
//                            sum[27] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{2},k{0},l{2};
//                            sum[28] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{2},k{0},l{1};
//                            sum[29] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }

//                        //sixth row
//                        {
//                            constexpr size_type i{0},j{1},k{0},l{0};
//                            sum[30] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{1},k{1},l{1};
//                            sum[31] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{1},k{2},l{2};
//                            sum[32] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{1},k{1},l{2};
//                            sum[33] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{1},k{0},l{2};
//                            sum[34] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                        {
//                            constexpr size_type i{0},j{1},k{0},l{1};
//                            sum[35] += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*weights;
//                        }
//                    }
//                }
//            }
//        }

//        data_(0,0) = sum[0]*eigth;
//        data_(0,1) = sum[1]*eigth;
//        data_(0,2) = sum[2]*eigth;
//        data_(0,3) = sum[3]*eigth;
//        data_(0,4) = sum[4]*eigth;
//        data_(0,5) = sum[5]*eigth;

//        data_(1,0) = sum[6]*eigth;
//        data_(1,1) = sum[7]*eigth;
//        data_(1,2) = sum[8]*eigth;
//        data_(1,3) = sum[9]*eigth;
//        data_(1,4) = sum[10]*eigth;
//        data_(1,5) = sum[11]*eigth;

//        data_(2,0) = sum[12]*eigth;
//        data_(2,1) = sum[13]*eigth;
//        data_(2,2) = sum[14]*eigth;
//        data_(2,3) = sum[15]*eigth;
//        data_(2,4) = sum[16]*eigth;
//        data_(2,5) = sum[17]*eigth;

//        data_(3,0) = sum[18]*2*eigth;
//        data_(3,1) = sum[19]*2*eigth;
//        data_(3,2) = sum[20]*2*eigth;
//        data_(3,3) = sum[21]*2*eigth;
//        data_(3,4) = sum[22]*2*eigth;
//        data_(3,5) = sum[23]*2*eigth;

//        data_(4,0) = sum[24]*2*eigth;
//        data_(4,1) = sum[25]*2*eigth;
//        data_(4,2) = sum[26]*2*eigth;
//        data_(4,3) = sum[27]*2*eigth;
//        data_(4,4) = sum[28]*2*eigth;
//        data_(4,5) = sum[29]*2*eigth;

//        data_(5,0) = sum[30]*2*eigth;
//        data_(5,1) = sum[31]*2*eigth;
//        data_(5,2) = sum[32]*2*eigth;
//        data_(5,3) = sum[33]*2*eigth;
//        data_(5,4) = sum[34]*2*eigth;
//        data_(5,5) = sum[35]*2*eigth;
//    }


//    template<typename Geometry, typename IntegrationPoints, typename Tensor>
//    constexpr inline auto eshelby_index_notation(size_type const i, size_type const j, size_type const k, size_type const l, IntegrationPoints const& points, Geometry const& a, Tensor const& L)const{
//        value_type sum{0};
//        uvw::math::matrix<value_type,3,3> K;
//        for(size_type M{0};M<points.number_outer;++M){
//            const auto Zeta3 = points.roots_outer[M];
//            const auto W_M = points.weights_outer[M];
//            for(size_type N{0};N<points.number_inner;++N){
//                const auto omega_q = M_PI*(points.roots_inner[N]+1);
//                const auto W_N = points.weights_inner[N];
//                const std::array<value_type,3> Xi{std::sqrt(1-Zeta3*Zeta3)*std::cos(omega_q)/a[0],std::sqrt(1-Zeta3*Zeta3)*std::sin(omega_q)/a[1],Zeta3/a[2]};//Xi(a,omega_q,Zeta3)};
//                compute_K(K,L,Xi);
//                value_type const D{compute_D(K)};
//                for(size_type m{0};m<3;++m){
//                    for(size_type n{0};n<3;++n){
//                        sum += L(m,n,k,l)*(compute_G(i,m,j,n,K,D,Xi) +  compute_G(j,m,i,n,K,D,Xi))*W_M*W_N;
//                    }
//                }
//            }
//        }
//        return sum/8;
//    }


//private:

//    template<typename Matrix, typename Tensor, typename Vector>
//    static constexpr inline auto compute_K(Matrix & K, Tensor const& L, Vector const& Xi){
//        for(std::size_t i{0};i<3;++i){
//            for(std::size_t k{0};k<3;++k){
//                value_type sum = 0;
//                for(std::size_t j{0};j<3;++j){
//                    for(std::size_t l{0};l<3;++l){
//                        sum += L(i,j,k,l)*Xi[j]*Xi[l];
//                    }
//                }
//                K(i,k) = sum;
//            }
//        }
//    }

//    template<typename Matrix>
//    static constexpr inline auto compute_D(Matrix const& K){
//        value_type sum{0};
//        for(std::size_t i{0};i<3;++i){
//            for(std::size_t j{0};j<3;++j){
//                for(std::size_t k{0};k<3;++k){
//                    sum += levi_civita_symbol(i+1,j+1,k+1)*K(i,0)*K(j,1)*K(k,2);
//                }
//            }
//        }
//        return sum;
//    }

//    template<typename Matrix, typename Vector>
//    static constexpr inline auto compute_G(size_type const i, size_type const j, size_type const k, size_type const l, Matrix const& K, value_type const D, Vector const& Xi){
//        return compute_Nij(i,j,K)*Xi[k]*Xi[l]/D;
//    }

//    template<typename Matrix>
//    static constexpr inline auto compute_Nij(size_type const i, size_type const j, Matrix const& K){
//        value_type N{0};
//        for(std::size_t k{0};k<3;++k){
//            for(std::size_t l{0};l<3;++l){
//                for(std::size_t m{0};m<3;++m){
//                    for(std::size_t n{0};n<3;++n){
//                        N += levi_civita_symbol(i+1,k+1,l+1)*levi_civita_symbol(j+1,m+1,n+1)*K(k,m)*K(l,n);
//                    }
//                }
//            }
//        }
//        return N*0.5;
//    }

//    static constexpr inline int levi_civita_symbol(size_type const i, size_type const j, size_type const k){
//        return -(int(j)-int(i))*(int(k)-int(j))*(int(i)-int(k))*0.5;
//    }

//    uvw::math::matrix<value_type,6,6> data_;
//};

#endif // ESHELBY_TENSOR_SOLID_ANISOTROPIC_BONES_H
