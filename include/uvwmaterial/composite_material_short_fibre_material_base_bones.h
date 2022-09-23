#ifndef COMPOSITE_MATERIAL_SHORT_FIBRE_MATERIAL_BASE_BONES_H
#define COMPOSITE_MATERIAL_SHORT_FIBRE_MATERIAL_BASE_BONES_H


template <typename _T, std::size_t _Dim>
class angular_central_gaussian_distribution
{
public:
    using value_type       = _T;
    using size_type        = std::size_t;
    using integration_data = SPHERE_LEBEDEV_RULE<value_type>;

    angular_central_gaussian_distribution(size_type __precision = 7):
        _N(),
        _B(),
        _data(),
        _weights(),
        _precision(__precision)
    {}

    constexpr inline auto & distribution_orientation_tensor(){
        return _N;
    }

    constexpr inline auto evaluate(){
        if(_precision == 0){
            throw std::runtime_error("angular_central_gaussian_distribution::evaluate(): precision == 0");
        }

        set_up_B();

        const auto size{integration_data::get_size(_precision)};
        _data.reserve(size/2ul + size*0.1);
        _weights.reserve(size/2ul + size*0.1);

        if constexpr (_Dim == 2){
            evaluate_2D();
        }else if constexpr (_Dim == 3){
            evaluate_3D();
        }

        //scale weights to get sum_i w_i = 1
        const auto sum_weights{std::accumulate(_weights.begin(), _weights.end(), value_type(0))};
        std::transform(_weights.begin(), _weights.end(), _weights.begin(),
                       [sum_weights](value_type const& value){return value/sum_weights;});

        _data.shrink_to_fit();
        _weights.shrink_to_fit();
    }

    constexpr inline auto get_directions_weights_data()const{
        return std::make_tuple(std::ref(_data), std::ref(_weights));
    }

    constexpr inline auto const& directions()const{
        return _data;
    }

    constexpr inline auto const& weights()const{
        return _weights;
    }

    constexpr inline auto number_of_integration_points()const{
        return _weights.size();
    }

    constexpr inline auto const& integration_precision()const{
        return _precision;
    }

    constexpr inline auto & integration_precision(){
        return _precision;
    }

private:
    constexpr inline auto evaluate_2D(){
        throw std::runtime_error("angular_central_gaussian_distribution::evaluate_2D(): not implemented jet");
    }

    inline auto evaluate_3D(){
        const auto size{integration_data::get_size(_precision)};
        const value_type* phi{integration_data::get_phi(_precision)};
        const value_type* theta{integration_data::get_theta(_precision)};
        const value_type* weights{integration_data::get_weight(_precision)};

        const value_type b1_sqrt{std::sqrt(_B[0])};
        const value_type b2_sqrt{std::sqrt(_B[1])};
        const value_type b3_sqrt{std::sqrt(_B[2])};

        for(std::size_t i{0}; i<size; ++i){
            //if (phi[i] >= 0 && phi[i] <= 90){
                value_type _x{std::cos((theta[i])*M_PI/180.)*std::sin(phi[i]*M_PI/180.)};
                value_type _y{std::sin((theta[i])*M_PI/180.)*std::sin(phi[i]*M_PI/180.)};
                value_type _z{std::cos(phi[i]*M_PI/180.)};
                _x = (std::abs(_x) < 1e-10 ? 0 : _x);
                _y = (std::abs(_y) < 1e-10 ? 0 : _y);
                _z = (std::abs(_z) < 1e-10 ? 0 : _z);
                const value_type norm{std::sqrt(std::pow(_x/b1_sqrt, 2) + std::pow(_y/b2_sqrt, 2) + std::pow(_z/b3_sqrt, 2))};
                _x /= (b1_sqrt*norm);
                _y /= (b2_sqrt*norm);
                _z /= (b3_sqrt*norm);
                if(_x >= 0){
                    _data.push_back(tmech::tensor<value_type, _Dim, 1>{_x, _y, _z});
                    _weights.push_back(weights[i]);
                }

            //}
        }
    }

    constexpr inline auto set_up_B(){
        if constexpr (_Dim == 2){
            set_up_B_2D();
        }else if constexpr (_Dim == 3){
            set_up_B_3D();
        }
    }

    constexpr inline auto set_up_B_2D(){
        throw std::runtime_error("angular_central_gaussian_distribution::set_up_B_2D(): not implemented jet");
    }

    constexpr inline auto set_up_B_3D(){
        //auto eigN{tmech::eigen_decomposition(_N)};
        //eigN.decompose();
        //_permut = eigN.permutation();

        //value_type a1{eigN.eigenvalues()[_permut[2]]};
        //value_type a2{eigN.eigenvalues()[_permut[1]]};
        //value_type a3{eigN.eigenvalues()[_permut[0]]};

        const value_type a1{_N(0,0)};
        const value_type a2{_N(1,1)};
        const value_type a3{_N(2,2)};

        if(a1 == 0 || a2 == 0 || a3 == 0){
            throw std::runtime_error("angular_central_gaussian_distribution::set_up_B_3D(): zero eigenvalue");
        }

        if(!almost_equal(a1+a2+a3, 1.0, 1)){
            throw std::runtime_error("angular_central_gaussian_distribution::set_up_B_3D(): trace(N)=" + std::to_string(a1+a2+a3) + " != 1" );
        }

        _B[0] = std::sqrt(a1);
        _B[1] = std::sqrt(a2);
        _B[2] = 1./_B[0]*_B[1];

        const value_type alpha{1};
        size_type method{0};
        for(; method<3; ++method){
            size_type iter{0}, max_iter{10000};
            value_type norm{0};
            for(; iter<max_iter; ++iter){
                const value_type a1_new{boost::math::ellint_rd(_B[1], _B[2], _B[0])/3.};
                const value_type a2_new{boost::math::ellint_rd(_B[0], _B[2], _B[1])/3.};
                const value_type a3_new{boost::math::ellint_rd(_B[0], _B[1], _B[2])/3.};

                const value_type R1{a1-a1_new};
                const value_type R2{a2-a2_new};
                const value_type R3{a3-a3_new};

                //std::cout<<"Eig A iter: "<<a1_new<<" "<<a2_new<<" "<<a3_new<<std::endl;
                switch (method) {
                case 0:
                {
                    _B[0] = std::abs(_B[0] - alpha*R1);
                    _B[1] = std::abs(_B[1] - alpha*R2);
                    _B[2] = 1./(_B[0]*_B[1]);
                    break;
                }
                case 1:
                {
                    _B[0] = std::abs(_B[0] - alpha*R1);
                    _B[2] = std::abs(_B[2] - alpha*R3);
                    _B[1] = 1./(_B[0]*_B[2]);
                    break;
                }
                case 2:
                {
                    _B[1] = std::abs(_B[1] - alpha*R2);
                    _B[2] = std::abs(_B[2] - alpha*R3);
                    _B[0] = 1./(_B[1]*_B[2]);
                    break;
                }
                }
                norm = std::sqrt(R1*R1 + R2*R2 + R3*R3);
                //std::cout<<"Number of iter "<<iter<<" norm "<<norm<<std::endl;
                if (norm <= 1e-14){
                    //std::cout<<"Number of iter "<<iter<<std::endl;
                    break;
                }
            }
            //std::cout<<b1<<" "<<b2<<" "<<b3<<std::endl;
            if (norm <= 1e-14){
                //std::cout<<"Number of iter "<<iter<<std::endl;
                break;
            }

            if(iter == max_iter-1 && method == 2){
                throw std::runtime_error("angular_central_gaussian_distribution::set_up_B_3D(): number of max iterations reached");
            }
        }
        std::cout<<_B[0]<<" "<<a1<<" "<<_B[1]<<" "<<a2<<" "<<_B[2]<<" "<<a3<<std::endl;
//        _B_tensor.fill(0);
//        for(size_type i{0}; i<3; ++i){
//            _B_tensor(i,i) = _B[i];
//        }
    }

protected:
    tmech::tensor<value_type, _Dim, 2> _N;
    std::array<value_type, _Dim> _B;
    tmech::tensor<value_type, _Dim, 2> _B_tensor;
    std::array<size_type, _Dim> _permut;
    std::vector<tmech::tensor<value_type, 3, 1>> _data;
    std::vector<value_type> _weights;
private:
    size_type _precision;

};




template <typename _T, std::size_t _Dim, typename _Container>
class short_fibre_composite_base
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    short_fibre_composite_base():
        _orientation_distribution()
    {}

    virtual ~short_fibre_composite_base(){}

    inline void init(){
        if constexpr (_Dim == 2){
            throw std::runtime_error("short_fibre_composite_base::init(): Dim == 2 not implemented jet");
        }else{
            _orientation_distribution.evaluate();
            _composites.reserve(_orientation_distribution.number_of_integration_points());
        }
    }

    constexpr inline auto number_of_composites()const{
        return _orientation_distribution.number_of_integration_points();
    }

    constexpr inline auto push_back(composite_material_base<_T, _Dim, _Container>* __composite){
        _composites.push_back(__composite);
    }

    template<typename _Derived>
    constexpr inline auto set_fibre_orientation_distribution_tensor(_Derived const& N){
        _orientation_distribution.distribution_orientation_tensor() = N;
    }

    constexpr inline auto const& directions()const{
        return _orientation_distribution.directions();
    }

    constexpr inline auto const& weights()const{
        return _orientation_distribution.weights();
    }

    constexpr inline auto const& composite_materials()const{
        return _composites;
    }

    constexpr inline auto& composite_materials(){
        return _composites;
    }

    constexpr inline auto const& integration_precision()const{
        return _orientation_distribution.integration_precision();
    }

    constexpr inline auto & integration_precision(){
        return _orientation_distribution.integration_precision();
    }

    virtual inline void update_strain() = 0;

protected:
    angular_central_gaussian_distribution<_T, _Dim> _orientation_distribution;
    std::vector<composite_material_base<_T, _Dim, _Container>*> _composites;
};

#endif // COMPOSITE_MATERIAL_SHORT_FIBRE_MATERIAL_BASE_BONES_H
