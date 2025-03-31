#ifndef UGPLUGIN_XBRAIDFORUG4_INITIALIZER_RANDOM_VALUE_INITIALIZER_HPP
#define UGPLUGIN_XBRAIDFORUG4_INITIALIZER_RANDOM_VALUE_INITIALIZER_HPP

#include "interface/initializer.hpp"
#include <random>




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class RandomValueInitializer : public BraidInitializer<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_VectorValueType =  typename TAlgebra::vector_type::value_type ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        //--------------------------------------------------------------------------------------------------------------
        int mode = 0;

        std::mt19937 generator;
        //std::default_random_engine generator;

        double param_0;
        double param_1;

        SP_DomainDisc domain_disc;
        //--------------------------------------------------------------------------------------------------------------

        RandomValueInitializer() {
            std::random_device rd; // seed random device
            this->generator = std::mt19937 (rd());
            param_0 = -1;
            param_1 =  1;
        };

        ~RandomValueInitializer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void initialize(SP_GridFunction& u, number time) override {
            if (time == this->t_start) {
                u = this->m_u0->clone();
            } else {
                u = this->m_u0->clone();
                if     ( mode == 0) {
                    this->fill_gridfunction_uniform(u);
                    //this->domain_disc->adjust_solution((*u),u->grid_level());
                }
                else if(mode == 1 ) { }
            }
        }

        void fill_gridfunction_uniform(SP_GridFunction& u_ref) {
            size_t szVector = u_ref->size();
            auto distribution = std::uniform_real_distribution<double>(param_0, param_1);

            for (size_t i = 0; i < szVector; i++) {
                u_ref->operator[](i) = distribution(generator);
            }
        }

        void fill_gridfunction_normal(SP_GridFunction& u_ref) {
            size_t szVector = u_ref->size();
            auto distribution = std::normal_distribution<double>(param_0, param_1);

            for (size_t i = 0; i < szVector; i++) {
                u_ref->operator[](i) = distribution(generator);
            }

        }
        //--------------------------------------------------------------------------------------------------------------
        void set_parameter_uniform(double min = -1.0, double max = 1.0) {
            this->mode = 0;
            this->param_0 = min;
            this->param_1 = max;
        }

        void set_parameter_normal(double mean = 0.0, double std = 1.0){
            this->mode = 1;
            this->param_0 = mean;
            this->param_1 = std;
        }

        void set_domain(SP_DomainDisc domain_disc) {
            this->domain_disc = domain_disc;
        }


        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif