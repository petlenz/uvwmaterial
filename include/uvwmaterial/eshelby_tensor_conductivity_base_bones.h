/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_CONDUCTIVITY_BASE_BONES_H
#define ESHELBY_TENSOR_CONDUCTIVITY_BASE_BONES_H

template <typename _T, std::size_t _Dim>
class eshelby_tensor_conductivity_base : public eshelby_tensor_base<_T, _Dim>
{
public:
    using value_type = _T;

    eshelby_tensor_conductivity_base();

    virtual ~eshelby_tensor_conductivity_base();

    inline virtual void init() = 0;

    inline virtual void reinit() = 0;

    inline constexpr auto const& tensor()const;

protected:
    tmech::tensor<value_type, _Dim, 2> _S;
};

#endif // ESHELBY_TENSOR_CONDUCTIVITY_BASE_BONES_H
