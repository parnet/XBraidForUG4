#ifndef UGPLUGIN_XBRAIDFORUG4_FACTORY_FIXED_STEP_THETA_INTEGRATOR_FACTORY_H
#define UGPLUGIN_XBRAIDFORUG4_FACTORY_FIXED_STEP_THETA_INTEGRATOR_FACTORY_H

#include "integrator/theta_conststep_integrator.hpp"
#include "interface/integrator_factory.hpp"




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class FixedStepThetaIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_TimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_ThetaIntegrator = ThetaConstStepIntegrator<TDomain, TAlgebra> ;
        using SP_ThetaIntegrator = SmartPtr<T_ThetaIntegrator> ;


        //--------------------------------------------------------------------------------------------------------------

        FixedStepThetaIntegratorFactory() = default;

        ~FixedStepThetaIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------
        void set_theta(double theta) {
            this->theta_ = theta;
        }

        void set_level_theta(int level, double theta) {
            this->level_theta_[level] = theta;
        }

        void set_num_steps(int steps) {
            this->num_steps_ = steps;
        }

        void set_level_num_steps(int level, int steps) {
            this->level_num_steps_[level] = steps;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(domain_disc_);
            integrator->set_solver(linear_solver_);
            integrator->set_theta(this->level_theta_[level]);
            integrator->set_num_steps(this->level_num_steps_[level]);
            return integrator;
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(domain_disc_);
            integrator->set_solver(linear_solver_);
            integrator->set_theta(this->theta_);
            return integrator;
        };

        void set_domain(SP_DomainDisc domain) {
            this->domain_disc_ = domain;
        }

        void set_solver(SP_Solver solver) {
            this->linear_solver_ = solver;
        }
        //--------------------------------------------------------------------------------------------------------------
        double theta_ = 1;
        int num_steps_ = 1;
        SP_Solver linear_solver_;
        SP_DomainDisc domain_disc_;
        std::vector<double> level_theta_;
        std::vector<double> level_num_steps_;
        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif