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

        RandomValueInitializer() {
            std::random_device rd; // seed random device
            this->generator_ = std::mt19937 (rd());
            param_0_ = -1;
            param_1_ =  1;
        };

        ~RandomValueInitializer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void initialize(SP_GridFunction& u, number time) override {
            if (time == this->t_start_) {
                u = this->u0_->clone();
            } else {
                u = this->u0_->clone();
                if     ( mode_ == 0) {
                    this->fill_gridfunction_uniform(u);
                    //this->domain_disc->adjust_solution((*u),u->grid_level());
                }
                else if(mode_ == 1 ) { }
            }
        }

        void fill_gridfunction_uniform(SP_GridFunction& u_ref) {
            size_t sz_vector = u_ref->size();
            auto distribution = std::uniform_real_distribution<double>(param_0_, param_1_);

            for (size_t i = 0; i < sz_vector; i++) {
                u_ref->operator[](i) = distribution(generator_);
            }
        }

        void fill_gridfunction_normal(SP_GridFunction& u_ref) {
            size_t sz_vector = u_ref->size();
            auto distribution = std::normal_distribution<double>(param_0_, param_1_);

            for (size_t i = 0; i < sz_vector; i++) {
                u_ref->operator[](i) = distribution(generator_);
            }

        }
        //--------------------------------------------------------------------------------------------------------------
        void set_parameter_uniform(double min = -1.0, double max = 1.0) {
            this->mode_ = 0;
            this->param_0_ = min;
            this->param_1_ = max;
        }

        void set_parameter_normal(double mean = 0.0, double std = 1.0){
            this->mode_ = 1;
            this->param_0_ = mean;
            this->param_1_ = std;
        }

        //2025-04 void set_domain(SP_DomainDisc domain_disc) {
        //2025-04     this->domain_disc_ = domain_disc;
        //2025-04 }

        //--------------------------------------------------------------------------------------------------------------
        int mode_ = 0;

        std::mt19937 generator_;
        //std::default_random_engine generator;

        double param_0_;
        double param_1_;

        //2025-04 SP_DomainDisc domain_disc_;

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif