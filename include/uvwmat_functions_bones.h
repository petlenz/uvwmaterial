#ifndef UVWMAT_FUNCTIONS_BONES_H
#define UVWMAT_FUNCTIONS_BONES_H

//in grad
template<typename T>
constexpr inline auto rotation_x1(T const alpha){
    constexpr auto factor{M_PI/180.};
    const auto alpha_radian{factor*alpha};
    const auto sn{std::sin(alpha_radian)};
    const auto cs{std::cos(alpha_radian)};

    return tmech::tensor<T, 3, 2>{1, 0,   0,
        0, cs, -sn,
                0, sn,  cs};
}

//in grad
template<typename T>
constexpr inline auto rotation_x2(T const alpha){
    constexpr auto factor{M_PI/180.};
    const auto alpha_radian{factor*alpha};
    const auto sn{std::sin(alpha_radian)};
    const auto cs{std::cos(alpha_radian)};

    return tmech::tensor<T, 3, 2>{cs, 0, sn,
                0,  1, 0,
                -sn, 0, cs};
}

template<std::size_t Dim = 3, typename T>
constexpr inline auto rotation_x3(T const alpha){
    constexpr auto factor{180./M_PI};
    const auto alpha_radian{factor*alpha};
    const auto sn{std::sin(alpha_radian)};
    const auto cs{std::cos(alpha_radian)};

    if constexpr (Dim == 2){
        return tmech::tensor<T, 2, 2>{cs, -sn,
                    sn,  cs};

    }else if constexpr (Dim == 3) {
        return tmech::tensor<T, 3, 2>{cs, -sn, 0,
                    sn,  cs, 0,
                    0,   0,  1};
    }else{
        static_assert (Dim == 2 || Dim == 3,"rotation_x1: is only valid for 2D and 3D");
    }
}

template<typename T>
constexpr inline auto rotation_x1_x2_x3(T const alpha1, T const alpha2, T const alpha3){
    constexpr T factor{M_PI/180.};
    const auto alpha1_radian{factor*alpha1};
    const auto alpha2_radian{factor*alpha2};
    const auto alpha3_radian{factor*alpha3};


    const auto sn1{std::sin(alpha1_radian)};
    const auto cs1{std::cos(alpha1_radian)};
    const auto sn2{std::sin(alpha2_radian)};
    const auto cs2{std::cos(alpha2_radian)};
    const auto sn3{std::sin(alpha3_radian)};
    const auto cs3{std::cos(alpha3_radian)};
    return tmech::tensor<T, 3, 2>{cs2*cs3,             -cs2*sn3,              sn2,
                sn1*sn2*cs3+cs1*sn3, -sn1*sn2*sn3+cs1*cs3, -sn1*cs2,
                -cs1*sn2*cs3+sn1*sn3,  cs1*sn2*sn3+sn1*cs3,  cs1*cs2};
}


template<typename T, typename Derived>
constexpr inline auto rotation_of_tensor(T const alpha1, T const alpha2, T const alpha3, tmech::tensor_base<Derived> const& data_base){
    using value_type = typename Derived::value_type;

    static_assert (Derived::rank() == 4,"only for fourth order tensors");
    static_assert (Derived::dimension() == 3,"only for 3 dimensional tensors");
    typename tmech::detail::result_expression_type<Derived>::result_type tensor{data_base.convert()};
    if constexpr(std::experimental::is_detected<tmech::detail::has_evaluate, typename std::remove_reference<Derived>::type>::value){
        tensor.evaluate();
    }

    //ist nicht allgemein gültig....
    const auto R{rotation_x1_x2_x3(alpha1, alpha2, alpha3)};
    tmech::tensor<value_type, 3, 4> temp;

    for(std::size_t i{0}; i<3; ++i){
        for(std::size_t j{0}; j<3; ++j){
            for(std::size_t k{0}; k<3; ++k){
                for(std::size_t l{0}; l<3; ++l){
                    T sum{0};
                    for(std::size_t m{0}; m<3; ++m){
                        for(std::size_t n{0}; n<3; ++n){
                            for(std::size_t o{0}; o<3; ++o){
                                for(std::size_t p{0}; p<3; ++p){
                                    sum += R(i,m)*R(j,n)*R(k,o)*R(l,p)*tensor(m,n,o,p);
                                }
                            }
                        }
                    }
                    temp(i,j,k,l) = sum;
                }
            }
        }
    }
    return temp;
}


template<typename T>
constexpr inline auto rotation_x3_x2_x1(T const alpha1, T const alpha2, T const alpha3){
    constexpr auto factor{M_PI/180.};
    const auto alpha1_radian{factor*alpha1};
    const auto alpha2_radian{factor*alpha2};
    const auto alpha3_radian{factor*alpha3};


    const auto sn1{std::sin(alpha1_radian)};
    const auto cs1{std::cos(alpha1_radian)};
    const auto sn2{std::sin(alpha2_radian)};
    const auto cs2{std::cos(alpha2_radian)};
    const auto sn3{std::sin(alpha3_radian)};
    const auto cs3{std::cos(alpha3_radian)};
    return tmech::tensor<T, 3, 2>{cs2*cs3, sn1*sn2*cs3-cs1*sn3, cs1*sn2*cs3+sn1*sn3,
                cs2*sn3, sn1*sn2*sn3+cs1*cs3, cs1*sn2*sn3-sn1*cs3,
                -sn2,    sn1*cs2,             cs1*cs2};
}

//https://de.wikipedia.org/wiki/Lam%C3%A9-Konstanten
template<typename T>
constexpr inline auto K_E_lambda(T const E, T const lambda){
    return ((E+3*lambda) + std::sqrt((E+3*lambda)*(E+3*lambda) -4*lambda*E)) / 6.;
}

template<typename T>
constexpr inline auto K_E_mu(T const E, T const mu){
    return (E*mu)/(3*(3*mu-E));
}

template<typename T>
constexpr inline auto K_E_nu(T const E, T const nu){
    return E/(3*(1-2*nu));
}

template<typename T>
constexpr inline auto K_lambda_mu(T const lambda, T const mu){
    return lambda + 2*mu/3;
}

template<typename T>
constexpr inline auto K_lambda_nu(T const lambda, T const nu){
    return (lambda*(1+nu))/(3*nu);
}

template<typename T>
constexpr inline auto K_mu_nu(T const mu, T const nu){
    return (2*mu*(1+nu))/(3*(1-2*nu));
}

template<typename T>
constexpr inline auto K_mu_M(T const mu, T const M){
    return M-4*mu/3;
}

template<typename T>
constexpr inline auto E_K_lambda(T const K, T const lambda){
    return (9*K*(K-lambda))/(3*K-lambda);
}

template<typename T>
constexpr inline auto E_K_mu(T const K, T const mu){
    return 9*K*mu/(3*K+mu);
}

template<typename T>
constexpr inline auto E_K_nu(T const K, T const nu){
    return 3*K*(1-2*nu);
}

template<typename T>
constexpr inline auto E_lambda_mu(T const lambda, T const mu){
    return (mu*(3*lambda+2*mu))/(lambda+mu);
}

template<typename T>
constexpr inline auto E_lambda_nu(T const lambda, T const nu){
    return (lambda*(1+nu)*(1-2*nu))/(nu);
}

template<typename T>
constexpr inline auto E_mu_nu(T const mu, T const nu){
    return 2*mu*(1+nu);
}

template<typename T>
constexpr inline auto E_mu_M(T const mu, T const M){
    return (mu*(3*M-4*mu))/(M-mu);
}

template<typename T>
constexpr inline auto lambda_K_E(T const K, T const E){
    return (3*K*(3*K-E))/(9*K-E);
}

template<typename T>
constexpr inline auto lambda_K_mu(T const K, T const mu){
    return K -2*mu/3;
}

template<typename T>
constexpr inline auto lambda_K_nu(T const K, T const nu){
    return 3*K*nu/(1+nu);
}

template<typename T>
constexpr inline auto lambda_E_mu(T const E, T const mu){
    return mu*(E-2*mu)/(3*mu-E);
}

template<typename T>
constexpr inline auto lambda_E_nu(T const E, T const nu){
    return E*nu/((1+nu)*(1-2*nu));
}

template<typename T>
constexpr inline auto lambda_mu_nu(T const mu, T const nu){
    2*mu*nu/(1-2*nu);
}

template<typename T>
constexpr inline auto lambda_mu_M(T const mu, T const M){
    return M-2*mu;
}

template<typename T>
constexpr inline auto mu_K_E(T const K, T const E){
    return 3*K*E/(9*K-E);
}

template<typename T>
constexpr inline auto mu_K_lambda(T const K, T const lambda){
    return 3*(K-lambda)*0.5;
}

template<typename T>
constexpr inline auto mu_K_nu(T const K, T const nu){
    return 3*K*(1-2*nu)/(2*(1+nu));
}

template<typename T>
constexpr inline auto mu_E_lambda(T const E, T const lambda){
    return (E-3*lambda)+std::sqrt((E-3*lambda)*(E-3*lambda) + 8*lambda*E)*0.25;
}

template<typename T>
constexpr inline auto mu_E_nu(T const E, T const nu){
    return E/(2*(1+nu));
}

template<typename T>
constexpr inline auto mu_lambda_nu(T const lambda, T const nu){
    return lambda*(1-2*nu)/(2*nu);
}

template<typename T>
constexpr inline auto nu_K_E(T const K, T const E){
    return (3*K-E)/(6*K);
}

template<typename T>
constexpr inline auto nu_K_lambda(T const K, T const lambda){
    return lambda/(3*K-lambda);
}

template<typename T>
constexpr inline auto nu_K_mu(T const K, T const mu){
    return (3*K-2*mu)/(2*(3*K+mu));
}

template<typename T>
constexpr inline auto nu_lambda_E(T const lambda, T const E){
    return -(E+lambda) + std::sqrt((E+lambda)*(E+lambda)+8*lambda*lambda)/(4*lambda);
}

template<typename T>
constexpr inline auto nu_E_mu(T const E, T const mu){
    return E/2*mu - 1.;
}

template<typename T>
constexpr inline auto nu_lambda_mu(T const lambda, T const mu){
    return lambda/(2*(lambda+mu));
}

template<typename T>
constexpr inline auto nu_mu_M(T const mu, T const M){
    return (M-2*mu)/(2*M-2*mu);
}


template<typename T, typename Tensor>
constexpr inline auto make_entry_minor_symmetric(std::size_t const i, std::size_t const j, std::size_t const k, std::size_t const l, T const entry, Tensor & C){
    C(i,j,k,l) = C(j,i,k,l) = C(i,j,l,k) = C(j,i,l,k) = entry;
    //C(i,j,k,l) = C(j,i,k,l) = C(i,j,l,k) = entry;
}

template<typename T, typename Tensor>
constexpr inline auto make_entry_major_symmetric(std::size_t const i, std::size_t const j, std::size_t const k, std::size_t const l, T const entry, Tensor & C){
    C(i,j,k,l) = C(k,l,i,j) = entry;
}

template<typename Derived>
constexpr inline auto extract_elasticity_parameter_K_mu(tmech::tensor_base<Derived> const& data_base){
    using value_type = typename Derived::value_type;
    typename tmech::detail::result_expression_type<Derived>::result_type tensor{data_base.convert()};
    if constexpr(std::experimental::is_detected<tmech::detail::has_evaluate, typename std::remove_reference<Derived>::type>::value){
        tensor.evaluate();
    }
    static_assert (Derived::rank() == 4,"extract_elasticity_parameter_K_mu: only for fourth-order tensors.");
    value_type K{0}, mu{0};
    for(std::size_t i{0}; i<Derived::dimension(); ++i){
        for(std::size_t j{0}; j<Derived::dimension(); ++j){
            K += tensor(i,i,j,j);
            mu += tensor(i,j,j,i);
        }
    }
    mu = (mu - K/3.)/10.;
    K = K/9.;
    return std::make_tuple(K, mu);
}

#endif // UVWMAT_FUNCTIONS_BONES_H
