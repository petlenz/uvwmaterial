#ifndef ESHELBY_TENSOR_CONDUCTIVITY_SPHERE_BONES_H
#define ESHELBY_TENSOR_CONDUCTIVITY_SPHERE_BONES_H

//Equivalent inclusion method for steady state heat conduction in composites
//Hiroshi Hatta Minoru Taya
//https://doi.org/10.1016/0020-7225(86)90011-X

//Random Heterogeneous Materials
//Microstructure and Macroscopic Properties
//Salvatore Torquato




template <typename _T, std::size_t _Dim>
class eshelby_tensor_conductivity_sphere :
        public eshelby_tensor_conductivity_base<_T, _Dim>
{
public:
    using value_type = _T;

    eshelby_tensor_conductivity_sphere(){}

    virtual ~eshelby_tensor_conductivity_sphere(){}

    virtual inline void init() override{
        if(!this->_is_init){
            eshelby_tensor_conductivity_sphere_func(this->_S);
            this->_is_init = true;
        }
    }

    inline virtual void reinit() {

    }
};


template <typename _T, std::size_t _Dim>
class eshelby_tensor_conductivity_cylindrical_ellipsoid :
        public eshelby_tensor_conductivity_base<_T, _Dim>
{
public:
    using value_type = _T;

    eshelby_tensor_conductivity_cylindrical_ellipsoid():
        _p()
    {}

    virtual ~eshelby_tensor_conductivity_cylindrical_ellipsoid(){}

    virtual inline void init() override{
        if(!this->_is_init){
            if(this->_parameter.size() != 2){
                throw std::runtime_error("eshelby_tensor_conductivity_general::init() number of parameters does not match");
            }
            const auto a{this->_parameter[0]};
            const auto b{this->_parameter[1]};
            if constexpr (_Dim == 3){
                eshelby_tensor_conductivity_cylindrical_ellipsoid_func(this->_S, _p, a, b);
            }else{
                tmech::tensor<value_type, 3, 2> Stemp;
                eshelby_tensor_conductivity_cylindrical_ellipsoid_func(Stemp, _p, a, b);
                this->_S = convert_3D_to_2D(Stemp);
            }
            this->_is_init = true;
        }
    }

    inline virtual void reinit() {

    }

    constexpr inline auto const& direction()const{
        return _p;
    }

    constexpr inline auto& direction(){
        return _p;
    }

protected:
    tmech::tensor<value_type, 3, 1> _p;
};


template <typename _T, std::size_t _Dim>
class eshelby_tensor_conductivity_cylinder :
        public eshelby_tensor_conductivity_base<_T, _Dim>
{
public:
    using value_type = _T;

    eshelby_tensor_conductivity_cylinder():
        _p()
    {}

    virtual ~eshelby_tensor_conductivity_cylinder(){}

    virtual inline void init() override{
        if(!this->_is_init){
            if constexpr (_Dim == 3){
                eshelby_tensor_conductivity_cylinder_func(this->_S, _p);
            }else{
                tmech::tensor<value_type, 3, 2> Stemp;
                eshelby_tensor_conductivity_cylinder_func(Stemp, _p);
                this->_S = convert_3D_to_2D(Stemp);
            }
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){

    }

    constexpr inline auto const& direction()const{
        return _p;
    }

    constexpr inline auto& direction(){
        return _p;
    }

protected:
    tmech::tensor<value_type, 3, 1> _p;
};


template <typename _T, std::size_t _Dim>
class eshelby_tensor_conductivity_general :
        public eshelby_tensor_conductivity_base<_T, _Dim>
{
public:
    using value_type = _T;

    eshelby_tensor_conductivity_general(){}

    virtual ~eshelby_tensor_conductivity_general(){}

    virtual inline void init() override{
        if(!this->_is_init){
            if(this->_parameter.size() != 3){
                throw std::runtime_error("eshelby_tensor_conductivity_general::init() number of parameters does not match");
            }
            const auto a{this->_parameter[0]};
            const auto b{this->_parameter[1]};
            const auto c{this->_parameter[2]};
            if constexpr (_Dim == 3){
                eshelby_tensor_conductivity_general_func(this->_S, a, b, c);
            }else{
                tmech::tensor<value_type, 3, 2> Stemp;
                eshelby_tensor_conductivity_general_func(Stemp, a, b, c);
                this->_S = convert_3D_to_2D(Stemp);
            }
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){
        this->_is_init = false;
        init();
    }
};
#endif // ESHELBY_TENSOR_CONDUCTIVITY_SPHERE_BONES_H
