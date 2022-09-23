#ifndef NONLOCAL_MATERIAL_BASE_BONES_H
#define NONLOCAL_MATERIAL_BASE_BONES_H

template <typename T, std::size_t Dim>
class nonlocal_material_base
{
public:
    using value_type = T;
    using size_type = std::size_t;
    using tensor2 = tmech::tensor<value_type, Dim, 2>;
    using key = size_type;

    nonlocal_material_base();

//    nonlocal_material_base(size_type const __number_of_nonlocal_variables);

//    nonlocal_material_base(std::vector<value_type> const& __nonlocal_variables);

    virtual ~nonlocal_material_base(){}

    virtual inline void update_nonlocal_variables() = 0;

    constexpr inline auto const& nonlocal_variables()const;

    constexpr inline auto& nonlocal_variables();

    template<typename _Key>
    constexpr inline auto number_of_variables(_Key const& __nonlocal_domain)const{
        //return _nonlocal_variables[__nonlocal_domain].size();
        return _nonlocal_variables.find(__nonlocal_domain)->second.size();
    }

    template<typename _Key>
    constexpr inline auto set_nonlocal_domain(_Key && __nonlocal_domain, size_type const __number_of_variables){
        _nonlocal_domains.push_back(__nonlocal_domain);
        _nonlocal_variables[__nonlocal_domain] = std::vector<value_type>(__number_of_variables);
        _source[__nonlocal_domain] = std::vector<tensor2>(__number_of_variables);
        _receiver[__nonlocal_domain] = std::vector<tensor2>(__number_of_variables);
        _element_marker[__nonlocal_domain] = false;
        _local_transformation[__nonlocal_domain] = tensor2();
        _internal_lengths[__nonlocal_domain] = std::vector<value_type>();
//        _nonlocal_domains.push_back(__nonlocal_domain);
//        _nonlocal_variables.insert({__nonlocal_domain, std::vector<value_type>(__number_of_variables)});
//        _source.insert({__nonlocal_domain, std::vector<tensor2>(__number_of_variables)});
//        _receiver.insert({__nonlocal_domain, std::vector<tensor2>(__number_of_variables)});
//        _element_marker.insert({__nonlocal_domain, false});
    }

    constexpr inline auto const& source_tensors()const{
        return _source;
    }

    constexpr inline auto const& receiver_tensors()const{
        return _receiver;
    }

    constexpr inline auto const& nonlocal_domains()const{
        return _nonlocal_domains;
    }

    template<typename _Key>
    constexpr inline auto const& source_tensor(_Key const& __nonlocal_domain)const{
        //assert(_source.find(__nonlocal_domain) != _source.end());
        return _source.find(__nonlocal_domain)->second;
        //return const_cast<boost::container::flat_map<key, std::vector<tmech::tensor<value_type, Dim, 2>>>&>(_source)[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto const& receiver_tensor(_Key const& __nonlocal_domain)const{
        //assert(_receiver.find(__nonlocal_domain) != _receiver.end());
        return _receiver.find(__nonlocal_domain)->second;
        //return const_cast<boost::container::flat_map<key, std::vector<tmech::tensor<value_type, Dim, 2>>>&>(_receiver)[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto& source_tensor(_Key const& __nonlocal_domain){
        //assert(_source.find(__nonlocal_domain) != _source.end());
        //return _source.find(__nonlocal_domain)->second;
        return _source[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto& receiver_tensor(_Key const& __nonlocal_domain){
        //assert(_receiver.find(__nonlocal_domain) != _receiver.end());
        //return _receiver.find(__nonlocal_domain)->second;
        return _receiver[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto const& nonlocal_variables(_Key const& __nonlocal_domain)const{
        return _nonlocal_variables.find(__nonlocal_domain)->second;
        //return const_cast<boost::container::flat_map<key, std::vector<value_type>>&>(_nonlocal_variables)[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto& nonlocal_variables(_Key const& __nonlocal_domain){
        //return _nonlocal_variables.find(__nonlocal_domain)->second;
        return _nonlocal_variables[__nonlocal_domain];
    }

    constexpr inline auto const& element_markers()const{
        return _element_marker;
    }

    template<typename _Key>
    constexpr inline auto const& element_marker(_Key const& __nonlocal_domain)const{
        return _element_marker.find(__nonlocal_domain)->second;
        //return const_cast<boost::container::flat_map<key, bool>&>(_element_marker)[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto& element_marker(_Key const& __nonlocal_domain){
        //return _element_marker.find(__nonlocal_domain)->second;
        return _element_marker[__nonlocal_domain];
    }

    constexpr inline auto const& internal_lengths()const{
        return _internal_lengths;
    }

    constexpr inline auto const& local_transformation()const{
        return _local_transformation;
    }

    template<typename _Key>
    constexpr inline auto const& internal_length(_Key const& __nonlocal_domain)const{
        return _internal_lengths.find(__nonlocal_domain)->second;
        //return const_cast<boost::container::flat_map<key, std::vector<value_type>>&>(_internal_lengths)[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto const& local_transformation(_Key const& __nonlocal_domain)const{
        return _local_transformation.find(__nonlocal_domain)->second;
        //return const_cast<boost::container::flat_map<key, tmech::tensor<value_type, Dim, 2>>&>(_local_transformation)[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto& internal_length(_Key const& __nonlocal_domain){
        return _internal_lengths[__nonlocal_domain];
    }

    template<typename _Key>
    constexpr inline auto& local_transformation(_Key const& __nonlocal_domain){
        return _local_transformation.find(__nonlocal_domain)->second;
    }
protected:
    boost::container::flat_map<key, bool> _element_marker;
    std::vector<key> _nonlocal_domains;
    boost::container::flat_map<key, std::vector<value_type>> _nonlocal_variables;
    boost::container::flat_map<key, std::vector<tmech::tensor<value_type, Dim, 2>>> _source;
    boost::container::flat_map<key, std::vector<tmech::tensor<value_type, Dim, 2>>> _receiver;
    boost::container::flat_map<key, std::vector<value_type>> _internal_lengths;
    std::map<key, tmech::tensor<value_type, Dim, 2>> _local_transformation;
//    std::unordered_map<key, std::vector<value_type>> _nonlocal_variables;
//    std::unordered_map<key, std::vector<tmech::tensor<value_type, Dim, 2>>> _source;
//    std::unordered_map<key, std::vector<tmech::tensor<value_type, Dim, 2>>> _receiver;
//    std::unordered_map<key, std::vector<value_type>> _internal_lengths;
//    std::unordered_map<key, tmech::tensor<value_type, Dim, 2>> _local_transformation;

    //std::unordered_map<std::string, bool> _element_marker;
    //std::vector<std::string> _nonlocal_domains;
//    std::unordered_map<std::string, std::vector<value_type>> _nonlocal_variables;
//    std::unordered_map<std::string, std::vector<tmech::tensor<value_type, Dim, 2>>> _source;
//    std::unordered_map<std::string, std::vector<tmech::tensor<value_type, Dim, 2>>> _receiver;
//    std::unordered_map<std::string, std::vector<value_type>> _internal_lengths;
//    std::unordered_map<std::string, tmech::tensor<value_type, Dim, 2>> _local_transformation;
};

#endif // NONLOCAL_MATERIAL_BASE_BONES_H
