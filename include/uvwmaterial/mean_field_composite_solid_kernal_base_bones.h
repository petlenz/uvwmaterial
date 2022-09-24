/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef MEAN_FIELD_COMPOSITE_SOLID_KERNAL_BASE_BONES_H
#define MEAN_FIELD_COMPOSITE_SOLID_KERNAL_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_solid_base;


template <typename _T, std::size_t _Dim, typename _Container>
class mean_field_composite_solid_kernal_base : public mean_field_composite_kernal_base<_T, _Dim, _Container>
{
public:
    mean_field_composite_solid_kernal_base() {}

    virtual ~mean_field_composite_solid_kernal_base(){}

    inline virtual void determine_strain_concentration_tensors(mean_field_composite_solid_base<_T, _Dim, _Container> & _composite) const = 0;
};

#endif // MEAN_FIELD_COMPOSITE_SOLID_KERNAL_BASE_BONES_H
