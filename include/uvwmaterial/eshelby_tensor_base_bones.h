/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef ESHELBY_TENSOR_BASE_BONES_H
#define ESHELBY_TENSOR_BASE_BONES_H

template <typename _T, std::size_t _Dim>
class eshelby_tensor_base
{
public:
    using value_type = _T;

    eshelby_tensor_base();

    virtual ~eshelby_tensor_base(){}

    inline virtual void init() = 0;

    inline virtual void reinit() = 0;

    constexpr inline bool is_initialized()const;

    template<typename ...Parameter>
    constexpr inline auto set_parameter(Parameter && ...__parameter);

    constexpr inline auto const& parameter()const{
        return _parameter;
    }

    constexpr inline auto& parameter(){
        return _parameter;
    }

protected:
    std::vector<value_type> _parameter;
    bool _is_init;
};

#endif // ESHELBY_TENSOR_BASE_BONES_H
