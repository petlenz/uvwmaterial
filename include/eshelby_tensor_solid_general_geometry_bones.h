#ifndef ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_BONES_H
#define ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_BONES_H

template <typename _T, std::size_t _Dim>
class eshelby_tensor_solid_general_geometry :
        public eshelby_tensor_solid_base<_T, _Dim>
{
    enum eshelby_type{SPHERE,
                      CYLINDER,
                      ELLIPTIC_CYLINDER,
                      SPHEROIDAL_1,
                      SPHEROIDAL_2,
                      SPHEROIDAL_3,
                      ELLIPSOID,
                      NON};
public:
    using value_type = _T;
    using size_type = std::size_t;

    eshelby_tensor_solid_general_geometry();

    virtual ~eshelby_tensor_solid_general_geometry();

    virtual inline void init() override;

    virtual inline void reinit()override{
        setup();
    }

private:
    constexpr inline auto setup();

    constexpr inline auto update_tensor();

    constexpr inline auto get_eshelby_type();

    constexpr inline auto sort();

    constexpr inline auto swap_elements();

    std::array<size_type, 3> _sort_index;
    eshelby_type _type;
};


#endif // ESHELBY_TENSOR_SOLID_GENERAL_GEOMETRY_BONES_H
