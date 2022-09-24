/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef TENSOR_ISOTROPIZATION_BASE_BONES_H
#define TENSOR_ISOTROPIZATION_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class tensor_isotropization_base
{
public:
    tensor_isotropization_base():
        _material_base(nullptr),
        _data(),
        _Iso(),
        _is_init(false)
    {}

    virtual ~tensor_isotropization_base(){}

    constexpr inline auto set_base_material(uvwmat::material_base<_T, _Dim, _Container>* __material_base){
        _material_base = __material_base;
    }

    inline virtual void init() = 0;

    inline virtual void update() = 0;

    inline virtual void update_tensor() = 0;

    constexpr inline auto const& data()const{
        return _data;
    }

    constexpr inline auto& data(){
        return _data;
    }

    constexpr inline auto const& tensor()const{
        return _Iso;
    }

    constexpr inline auto& tensor(){
        return _Iso;
    }

protected:
    material_base<_T, _Dim, _Container>* _material_base;
    std::vector<_T> _data;
    tmech::tensor<_T, _Dim, 4> _Iso;
    bool _is_init;
};

#endif // TENSOR_ISOTROPIZATION_BASE_BONES_H
