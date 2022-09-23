#ifndef COMPOSITE_MATERIAL_BASE_BONES_H
#define COMPOSITE_MATERIAL_BASE_BONES_H

template <typename _T, std::size_t _Dim, typename _Container>
class composite_material_base
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    composite_material_base();

    virtual ~composite_material_base();

    inline virtual value_type& volume_fraction(size_type const __idx) = 0;

    inline virtual value_type const& volume_fraction(size_type const __idx) const = 0;

    inline virtual size_type number_of_materials() const = 0;

    inline virtual material_base<_T, _Dim, _Container> const* material(size_type const __idx) const = 0;

    inline virtual material_base<_T, _Dim, _Container>* material(size_type const __idx) = 0;

    constexpr inline void init_materials(){
        _materials.resize(number_of_materials());
        for(size_type i{0}; i<_materials.size(); ++i){
            _materials[i] = material(i);
        }
    }

    constexpr inline auto const& materials()const{
        return _materials;
    }

    constexpr inline auto& materials(){
        return _materials;
    }



protected:
    bool _is_init;
    std::vector<material_base<_T, _Dim, _Container>*> _materials;
};



#endif // COMPOSITE_MATERIAL_BASE_BONES_H
