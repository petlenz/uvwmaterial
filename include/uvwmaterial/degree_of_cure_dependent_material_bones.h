#ifndef DEGREE_OF_CURE_DEPENDENT_MATERIAL_BONES_H
#define DEGREE_OF_CURE_DEPENDENT_MATERIAL_BONES_H

template <typename T>
class degree_of_cure_dependent_material_base
{
public:
    using value_type = T;

    degree_of_cure_dependent_material_base():
        _curing(0),
        _dcuring(0)
    {}

    virtual ~degree_of_cure_dependent_material_base(){}


    constexpr inline auto& degree_of_cure(){
        return _curing;
    }

    constexpr inline auto const& degree_of_cure()const{
        return _curing;
    }

    constexpr inline auto& delta_degree_of_cure(){
        return _dcuring;
    }

    constexpr inline auto const& delta_degree_of_cure()const{
        return _dcuring;
    }

protected:
    value_type _curing;
    value_type _dcuring;
};

#endif // DEGREE_OF_CURE_DEPENDENT_MATERIAL_BONES_H
