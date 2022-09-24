/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_MEAT_H
#define ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_MEAT_H

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_general_geometry<_T, _Dim>::eshelby_tensor_solid_general_geometry():
    _sort_index({0,1,2}),
    _type(eshelby_type::NON)
{}

template <typename _T, std::size_t _Dim>
eshelby_tensor_solid_general_geometry<_T, _Dim>::~eshelby_tensor_solid_general_geometry(){}

template <typename _T, std::size_t _Dim>
inline void eshelby_tensor_solid_general_geometry<_T, _Dim>::init(){
    if(!this->_is_init){
        setup();
        this->_is_init = true;
    }
}

template <typename _T, std::size_t _Dim>
constexpr inline auto eshelby_tensor_solid_general_geometry<_T, _Dim>::setup(){
    sort();
    get_eshelby_type();
    update_tensor();
    if(_sort_index != std::array<std::size_t, 3>{0,1,2}){//_sort_index[0] != 0 || _sort_index[1] != 1 || _sort_index[2] != 2){
        //std::cout<<"Swap elements "<<_sort_index[0]<<" "<<_sort_index[1]<<" "<<_sort_index[2]<<std::endl;
        swap_elements();
    }
    //value_type* _a{&this->_parameter[1]};
    //std::cout<<_a[_sort_index[0]]<<" "<<_a[_sort_index[1]]<<" "<<_a[_sort_index[2]]<<std::endl;
    if(false){
    value_type S_v[36];
    tmech::convert_tensor_to_voigt(this->_S, S_v);
    for(int i{0}; i<6; ++i){
        for(int j{0}; j<6; ++j){
            std::cout<<S_v[i*6+j]<<" ";
        }
        std::cout<<std::endl;
    }
    }
}

template <typename _T, std::size_t _Dim>
constexpr inline auto eshelby_tensor_solid_general_geometry<_T, _Dim>::update_tensor(){
    value_type* _a{&this->_parameter[1]};
    //_a[0] == a; _a[1] == b; _a[2] == c
    //a1 < a2 < a3
    const auto a1{_a[_sort_index[0]]}, a2{_a[_sort_index[1]]}, a3{_a[_sort_index[2]]};
    this->_S.fill(0);
    switch (_type) {
    case eshelby_type::SPHERE:
        eshelby_tensor_sphere_func(this->_S, this->_parameter[0]);
        break;
    case eshelby_type::CYLINDER:
        eshelby_tensor_elliptic_cylinder_func(this->_S, this->_parameter[0], a1, a2);
        break;
    case eshelby_type::ELLIPTIC_CYLINDER:
        eshelby_tensor_elliptic_cylinder_func(this->_S, this->_parameter[0], a1, a2);
        break;
    case eshelby_type::SPHEROIDAL_1://aligned x1 direction equal_23==true
        if constexpr (_Dim == 3){
            eshelby_tensor_spheroidal_func(this->_S, this->_parameter[0], a1/a2, tmech::tensor<_T,3,1>{1,0,0});
        }
        break;
    case eshelby_type::SPHEROIDAL_2://aligned x2 direction equal_13==true
        if constexpr (_Dim == 3){
            eshelby_tensor_spheroidal_func(this->_S, this->_parameter[0], a2/a1, tmech::tensor<_T,3,1>{0,1,0});
        }
        break;
    case eshelby_type::SPHEROIDAL_3://aligned x3 direction equal_12==true
        if constexpr (_Dim == 3){
            eshelby_tensor_spheroidal_func(this->_S, this->_parameter[0], a3/a2, tmech::tensor<_T,3,1>{0,0,1});
        }
        break;
    case eshelby_type::ELLIPSOID:
        //a1 < a2 < a3
        //aligned in x1 direction
        eshelby_tensor_ellipsoid_func(this->_S, this->_parameter[0], a3, a2, a1);
        break;
    case eshelby_type::NON:
        throw std::runtime_error("eshelby_tensor_general::update() no matching geometry type");
        break;
    }
}

template <typename _T, std::size_t _Dim>
constexpr inline auto eshelby_tensor_solid_general_geometry<_T, _Dim>::get_eshelby_type(){
    value_type* _a{&this->_parameter[1]};
    const auto a1{_a[_sort_index[0]]}, a2{_a[_sort_index[1]]}, a3{_a[_sort_index[2]]};
    const auto max_element{(a3 == std::numeric_limits<value_type>::infinity() ? std::max(a1,a2) : a3)};
    const auto eps{std::numeric_limits<value_type>::epsilon()*max_element};
    const auto equal_12{std::fabs(a1-a2) <= eps};
    const auto equal_13{std::fabs(a1-a3) <= eps};
    const auto equal_23{std::fabs(a2-a3) <= eps};
    //

    if(equal_12 && equal_13 && equal_23){
        //no need for to swapp element
        _type = eshelby_type::SPHERE;
        //std::cout<<"SPHERE "<<a1<<" "<<a2<<" "<<a3<<std::endl;
        return ;
    }

    //check if one pair is equal
    //this results into spheroidal inclusion with aspect ratio Ar = length/diameter
    //the spheroidal inclusion is algined in x1 direction
    //so the case equal_23 has the correct entries, for the other cases
    //entries are needed to be swapped
    if(equal_12){
        //aligned in x3 direction
        _type = eshelby_type::SPHEROIDAL_3;
        return;
    }

    if(equal_13){
        //aligned in x2 direction
        _type = eshelby_type::SPHEROIDAL_2;
        return;
    }

    if(equal_23){
        //aligned in x1 direction
        _type = eshelby_type::SPHEROIDAL_1;
        return;
    }

    if(a3 == std::numeric_limits<value_type>::infinity() || a3/a2 >= 10){
        if(equal_12){
            //std::cout<<"CYLINDER "<<a1<<" "<<a2<<" "<<a3<<std::endl;
            _type = eshelby_type::ELLIPTIC_CYLINDER;
            return ;
        }else{
            //std::cout<<"ELLIPTIC CYLINDER "<<a1<<" "<<a2<<" "<<a3<<std::endl;
            _type = eshelby_type::ELLIPTIC_CYLINDER;
            return ;
        }
    }

    //std::cout<<"General ellipsoid "<<a1<<" "<<a2<<" "<<a3<<std::endl;
    _type = eshelby_type::ELLIPSOID;
    return ;
}

template <typename _T, std::size_t _Dim>
constexpr inline auto eshelby_tensor_solid_general_geometry<_T, _Dim>::sort(){
    value_type* _a{&this->_parameter[1]};
    _sort_index = {0,1,2};
    if (_a[_sort_index[0]] > _a[_sort_index[1]])
        std::swap(_sort_index[0], _sort_index[1]);
    if (_a[_sort_index[0]] > _a[_sort_index[2]])
        std::swap(_sort_index[0], _sort_index[2]);
    if (_a[_sort_index[1]] > _a[_sort_index[2]])
        std::swap(_sort_index[1], _sort_index[2]);
}

template <typename _T, std::size_t _Dim>
constexpr inline auto eshelby_tensor_solid_general_geometry<_T, _Dim>::swap_elements(){
    const auto copyS{this->_S};
    auto& S{this->_S};
    for(std::size_t i{0}; i<_Dim; ++i){
        const auto I{_sort_index[i]};
        for(std::size_t j{0}; j<_Dim; ++j){
            const auto J{_sort_index[j]};
            for(std::size_t k{0}; k<_Dim; ++k){
                const auto K{_sort_index[k]};
                for(std::size_t l{0}; l<_Dim; ++l){
                    const auto L{_sort_index[l]};
                    S(i,j,k,l) = copyS(I,J,K,L);
                    //std::swap(S(I,J,K,L), S(i,j,k,l));
                }
            }
        }
    }
}
#endif // ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_MEAT_H
