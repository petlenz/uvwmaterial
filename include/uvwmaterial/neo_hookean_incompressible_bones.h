/***************************************************************************
* Copyright (c) Peter Lenz                                                 *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/
#ifndef NEO_HOOKEAN_INCOMPRESSIBLE_BONES_H
#define NEO_HOOKEAN_INCOMPRESSIBLE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class neo_hookean_incompressible :
        public finite_strain_solid_material_base<_T, _Dim>,
        public solid_material_base<_T, _Dim, _Container>
{
public:
    neo_hookean_incompressible() {}

    inline virtual void init(){
        if(!this->_is_init){
            this->_C.fill(0);
            this->_is_init = true;
        }
    }

    inline virtual void reinit(){

    }

    inline virtual void update(){

    }

    inline virtual void update_stress(){

    }

    inline virtual void update_tangent(){

    }


};
#endif // NEO_HOOKEAN_INCOMPRESSIBLE_BONES_H
