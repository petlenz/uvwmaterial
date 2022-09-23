#ifndef ESHELBY_TENSOR_FUNCTIONS_H
#define ESHELBY_TENSOR_FUNCTIONS_H



template<typename _Eshelby, typename _T>
constexpr inline auto eshelby_tensor_conductivity_general_func(_Eshelby & _S, _T const& _a1, _T const& _a2, _T const& _a3){
    const auto fac1{_a1*_a2*_a3/_T(2)};
    const auto aa1{_a1*_a1}, aa2{_a2*_a2}, aa3{_a3*_a3};
    _S(0,0) = fac1*_T(2)/_T(3)*boost::math::ellint_rj(aa1,aa2,aa3,aa1);
    _S(1,1) = fac1*_T(2)/_T(3)*boost::math::ellint_rj(aa1,aa2,aa3,aa2);
    _S(2,2) = fac1*_T(2)/_T(3)*boost::math::ellint_rj(aa1,aa2,aa3,aa3);
}

template<typename _Tensor_S>
static constexpr inline auto eshelby_tensor_conductivity_sphere_func(_Tensor_S & __S){
    using value_type = typename _Tensor_S::value_type;
    constexpr auto Dim{_Tensor_S::dimension()};
    __S = tmech::eye<value_type, Dim, 2>()/Dim;
}

template<typename _Tensor_S, typename _Direction>
static constexpr inline auto eshelby_tensor_conductivity_cylinder_func(_Tensor_S & __S, _Direction const& __p){
    using value_type = typename _Tensor_S::value_type;
    constexpr auto Dim{_Tensor_S::dimension()};
    const tmech::eye<value_type, Dim, 2> I;
    __S = 0.5*(I - tmech::otimes(__p, __p));
}

template<typename _Tensor_S, typename _Direction, typename _T>
static constexpr inline auto eshelby_tensor_conductivity_cylindrical_ellipsoid_func(_Tensor_S & __S, _Direction const& __p, _T const __a, _T const __b){
    using value_type = typename _Tensor_S::value_type;
    constexpr auto Dim{_Tensor_S::dimension()};
    const tmech::eye<value_type, Dim, 2> I;
    auto Q{[](value_type const __a, value_type const __b){
            const auto val1{__a*__a/(__b*__b)};
            if(__b/__a >= 1.0){
                const auto gamma_b{std::sqrt(1 - val1)};
                return 0.5*(1 + (1-std::log((1+gamma_b)/(1-gamma_b))/(2*gamma_b))/(val1-1));
            }else{
                const auto gamma_a{std::sqrt(val1 - 1)};
                return 0.5*(1 + (1-std::atan(gamma_a)/gamma_a)/(val1-1));
            }
        }};
    __S = Q(__a,__b)*(I - 2*tmech::otimes(__p, __p));
}

template<typename _Eshelby, typename _T>
constexpr inline auto eshelby_tensor_sphere_func(_Eshelby & _S, _T const& _nue){
    using value_type = typename _Eshelby::value_type;

    const tmech::eye<value_type, _Eshelby::dimension(), 2> I;
    const auto IIsym = 0.5*(tmech::otimesu(I,I)+tmech::otimesl(I,I));
    const auto II = tmech::otimes(I,I);
    const auto IIvol = II/3.0;

    const double gamma{(1+_nue)/(9*(1-_nue))};
    const double delta{(4-5*_nue)/(15*(1-_nue))};
    _S = gamma*II + delta*2*(IIsym-IIvol);
}

template<typename _Eshelby, typename _T>
constexpr inline auto eshelby_tensor_spheroids_func(_Eshelby & _S, _T const& _nue, _T const& _Ar){
    using value_type = typename _Eshelby::value_type;
    //https://upcommons.upc.edu/bitstream/handle/2117/22253/MONOGRAFIA_ORTOLANO_hDEZ.pdf?sequence=1
}

template<typename _Eshelby, typename _T>
constexpr inline auto eshelby_tensor_elliptic_cylinder_func(_Eshelby & _S, _T const& _nue, _T const& _a1, _T const& _a2){
    //a3 --> infty; a1/=a2
    //--> aligned in x3 diretion
    const auto val1{1/(2*(1-_nue))};
    const auto val2{(_a1+_a2)*(_a1+_a2)};
    const auto aa1{_a1*_a1};
    const auto aa2{_a2*_a2};

    _S(0,0,0,0) = val1*((aa2+2*_a1*_a2)/val2 + (1-2*_nue)*_a2/(_a1+_a2));
    _S(1,1,1,1) = val1*((aa1+2*_a1*_a2)/val2 + (1-2*_nue)*_a1/(_a1+_a2));
    _S(0,0,1,1) = val1*(aa2/val2 - (1-2*_nue)*_a2/(_a1+_a2));
    _S(1,1,0,0) = val1*(aa1/val2 - (1-2*_nue)*_a1/(_a1+_a2));
    _S(1,1,2,2) = val1*2*_nue*_a1/(_a1+_a2);
    _S(0,0,2,2) = val1*2*_nue*_a2/(_a1+_a2);
    _S(0,1,0,1) = val1*((aa1+aa2)/(2*val2) + (1-2*_nue)/2);
    _S(1,2,1,2) = _a1/(2*(_a1+_a2));
    _S(2,0,2,0) = _a2/(2*(_a1+_a2));
}


template<typename _Eshelby, typename _T>
constexpr inline auto eshelby_tensor_ellipsoid_func(_Eshelby & _S, _T const& _nue, _T const& _a1, _T const& _a2, _T const& _a3){
    //a1 > a2 > a3
    //aligned in x1 direction
    const auto theta{std::asin(std::sqrt(1.-_a3*_a3/(_a1*_a1)))};//(_a1*_a1 - _a3*_a3)/(_a1*_a1)
    const auto k{std::sqrt((_a1*_a1 - _a2*_a2)/(_a1*_a1 - _a3*_a3))};
    const auto F{boost::math::ellint_1(k, theta)};
    const auto E{boost::math::ellint_2(k, theta)};
    std::array<_T,3> I;
       I[0] = 4*M_PI*_a1*_a2*_a3*(F-E)/((_a1*_a1 - _a2*_a2)*std::sqrt((_a1*_a1 - _a3*_a3)));
       I[2] = 4*M_PI*_a1*_a2*_a3*((_a2*std::sqrt(_a1*_a1 - _a3*_a3)/(_a1*_a3)) - E)/((_a2*_a2 - _a3*_a3)*std::sqrt(_a1*_a1 - _a3*_a3));
       I[1] = 4*M_PI - (I[0]+I[2]);
    const std::array<_T,3> a{_a1,_a2,_a3};
    std::array<std::array<_T,3>,3> II{0};
    for(std::size_t i{0}; i<3; ++i){
        for(std::size_t j{0}; j<3; ++j){
            if(i!=j){
                II[i][j] = (I[j] - I[i])/(a[i]*a[i] - a[j]*a[j]);
            }
        }
    }

    for(std::size_t i{0}; i<3; ++i){
        _S(i,i,i,i) = (3*a[i]*a[i]*II[i][i] + (1-2*_nue)*I[i])/(8*M_PI*(1-_nue));
        for(std::size_t j{0}; j<3; ++j){
            if(j != i){
                _S(i,i,j,j) = (a[j]*a[j]*II[i][j] - (1-2*_nue)*I[i])/(8*M_PI*(1-_nue));
                _S(i,j,i,j) = ((a[i]*a[i] + a[j]*a[j])*II[i][j] + (1-2*_nue)*(I[i]+I[j]))/(16*M_PI*(1-_nue));
                _S(j,i,j,i) = _S(i,j,i,j);
                _S(i,j,j,i) = _S(i,j,i,j);
                _S(j,i,i,j) = _S(i,j,i,j);
            }
        }
    }
}


template<typename _Eshelby, typename _T, typename _Direction>
constexpr inline auto eshelby_tensor_spheroidal_func(_Eshelby & _S, _T const __nue, _T const __eta, _Direction const& __n){
    const _T eta2{__eta*__eta};
    //eta = length/diameter
    auto Tg{[__eta, eta2](){
            if(__eta < 1.0){
                return (__eta/(std::pow((1.0-eta2),1.5)))*(std::acos(__eta) - __eta*std::sqrt(1.0-eta2));
            }else{
                return (__eta/(std::pow(eta2-1.0, 1.5)))*(__eta*std::sqrt(eta2 - 1.0) - (std::acosh(__eta)));
            }
        }};

    const _T g{Tg()};

    const tmech::eye<_T, 3, 2> I;
    const auto p{tmech::otimes(__n,__n)};
    const auto II1{0.5*(tmech::otimesu(I,I) + tmech::otimesl(I,I))};
    const auto II2{tmech::otimes(I,I)};
    const auto II3{tmech::otimes(p,I)};
    const auto II4{tmech::otimes(I,p)};
    const auto II5{(tmech::otimesu(p,I) + tmech::otimesl(p,I) + tmech::otimesu(I,p) + tmech::otimesl(I,p))*0.5};
    const auto II6{tmech::otimes(p, p)};

    const _T c{(eta2 - 1)*__nue - eta2 + 1};
    const _T c1{((8*(eta2-1)*__nue - 4*eta2 + 7)*g - 2*eta2)/(8*c)};
    const _T c2{-((8*(eta2-1)*__nue - 4*eta2 + 1)*g + 2*eta2)/(16*c)};
    const _T c3{((24*(eta2-1)*__nue - 12*eta2 - 3)*g + 16*(1-eta2)*__nue + 10*eta2)/(16*c)};
    const _T c4{-((12*eta2+3)*g-10*eta2)/(16*c)};
    const _T c5{-((12*(eta2 - 1)*__nue + 15)*g+8*(1-eta2)*__nue - 2*eta2 - 8)/(8*c)};
    const _T c6{((60*eta2+45)*g - 54*eta2 - 16)/(16*c)};

    _S = c1*II1 + c2*II2 + c3*II3 + c4*II4 + c5*II5 + c6*II6;
}

#endif // ESHELBY_TENSOR_FUNCTIONS_H
