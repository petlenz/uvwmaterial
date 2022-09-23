#ifndef TENSOR_ISOTROPIZATION_GENERAL_BONES_H
#define TENSOR_ISOTROPIZATION_GENERAL_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class tensor_isotropization_general :
        public tensor_isotropization_base<_T, _Dim, _Container>
{
public:
    using value_type = _T;

    tensor_isotropization_general()
    {}

    virtual ~tensor_isotropization_general(){}

    inline virtual void init() override {
        if(!this->_is_init){
            if(!this->_material_base){
                throw std::runtime_error("eshelby_tensor_general_isotropization_base::init() base material is not set");
            }
            this->_parameter.resize(2);
            const tmech::eye<value_type,_Dim,2> I;
            //fourth order symmetric idenity tensor
            const auto IIsym = 0.5*(tmech::otimesu(I,I)+tmech::otimesl(I,I));
            //fourth order identity tensor
            const auto II = tmech::otimes(I,I);
            const auto IIvol = II/_Dim;
            this->_data[0] = update(IIvol)/_Dim;
            this->_is_init = true;
        }
    }

    inline virtual void update() override {
        const tmech::eye<value_type,_Dim,2> I;
        //fourth order symmetric idenity tensor
        const auto IIsym = 0.5*(tmech::otimesu(I,I)+tmech::otimesl(I,I));
        //fourth order identity tensor
        const auto II = tmech::otimes(I,I);
        const auto IIvol = II/_Dim;
        const auto IIdev = IIsym-IIvol;
        this->_data[1] = update(IIdev)/10.0;
    }

    inline virtual void update_tensor() override {
        const tmech::eye<value_type,_Dim,2> I;
        //fourth order symmetric idenity tensor
        const auto IIsym = 0.5*(tmech::otimesu(I,I)+tmech::otimesl(I,I));
        //fourth order identity tensor
        const auto II = tmech::otimes(I,I);
        const auto IIvol = II/3.0;
        const auto IIdev = IIsym-IIvol;
        this->_Iso = 3*this->_data[0]*IIvol + 2*this->_data[1]*IIdev;
    }

private:
    template<typename _Tensor>
    constexpr inline auto update(_Tensor const& _tensor){
        auto matrix = static_cast<uvwmat::solid_material_base<_T, _Dim, _Container>*>(this->_material_base);
        return tmech::ddcontract(_tensor, matrix.get()->tangent_tensor());
    }
};

#endif // TENSOR_ISOTROPIZATION_GENERAL_BONES_H
