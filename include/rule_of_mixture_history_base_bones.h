#ifndef RULE_OF_MIXTURE_HISTORY_BASE_BONES_H
#define RULE_OF_MIXTURE_HISTORY_BASE_BONES_H

template<typename _T, std::size_t _Dim, typename _Container>
class rule_of_mixture_history_base :
        public rule_of_mixture_composite_base<_T, _Dim, _Container>,
        public history_material_base<_T>
{
public:
    using value_type = _T;
    using size_type  = std::size_t;

    constexpr rule_of_mixture_history_base(){}

    virtual ~rule_of_mixture_history_base(){}

    constexpr inline void init(){
        rule_of_mixture_composite_base<_T, _Dim, _Container>::init();
        auto number_of_his_var{0};
        for(auto material : this->_materials){
            auto his_material{dynamic_cast<history_material_base<value_type>*>(material)};
            if(his_material){
                number_of_his_var += his_material->get_number_of_variables();
            }
        }
        this->history.resize(number_of_his_var);
    }

protected:
    constexpr inline auto set_local_history(){
        size_type iter{0};

        for(auto material : this->_materials){
            auto his_material{dynamic_cast<history_material_base<value_type>*>(material)};
            if(his_material){
                const auto his_number{his_material->get_number_of_variables()};
                his_material->set_history(&this->history[iter], &this->history[iter+his_number]);
                iter += his_number;
            }
        }
    }

    constexpr inline auto get_local_history(){
        std::vector<value_type> his;
        his.reserve(this->history.size());

        for(auto material : this->_materials){
            auto his_material{dynamic_cast<history_material_base<value_type>*>(material)};
            if(his_material){

                for(const auto& h : his_material->get_history()){
                    his.push_back(h);
                }
            }
        }

        this->history = his;
    }
};

#endif // RULE_OF_MIXTURE_HISTORY_BASE_BONES_H
