/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_SOLID_CYLINDRICAL_ELLIPSOID_BONES_H
#define ESHELBY_TENSOR_SOLID_CYLINDRICAL_ELLIPSOID_BONES_H

//https://sci-hub.mksa.top/10.1002/pc.750050413
//
template <typename _T, std::size_t _Dim>
class eshelby_tensor_solid_cylindrical_ellipsoid :
        public eshelby_tensor_solid_base<_T, _Dim>
{
public:
    using value_type = _T;
    using base_class = eshelby_tensor_solid_base<_T, _Dim>;

    eshelby_tensor_solid_cylindrical_ellipsoid();

    virtual ~eshelby_tensor_solid_cylindrical_ellipsoid();

    inline void init() override;

    inline virtual void reinit(){}

    constexpr inline auto& direction();

    constexpr inline auto const& direction()const;

private:
    template<typename _N>
    static constexpr inline auto ellipsoid_direction(_T const __nue, _T const __eta, _N const& __n);

    tmech::tensor<value_type, 3, 1> _n;
};


#endif // ESHELBY_TENSOR_SOLID_CYLINDRICAL_ELLIPSOID_BONES_H
