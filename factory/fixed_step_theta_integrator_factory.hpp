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
        double theta = 1;
        int num_steps = 1;
        SP_Solver m_linear_solver;
        SP_DomainDisc m_domain_disc;
        std::vector<double> level_theta;
        std::vector<double> level_num_steps;
        //--------------------------------------------------------------------------------------------------------------

        FixedStepThetaIntegratorFactory() = default;

        ~FixedStepThetaIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------
        void set_theta(double p_theta) {
            this->theta = p_theta;
        }

        void set_level_theta(int level, double p_theta) {
            this->level_theta[level] = p_theta;
        }

        void set_num_steps(int steps) {
            this->num_steps = steps;
        }

        void set_level_num_steps(int level, int steps) {
            this->level_num_steps[level] = steps;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_theta(this->level_theta[level]);
            integrator->set_num_steps(this->level_num_steps[level]);
            return integrator;
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_theta(this->theta);
            return integrator;
        };

        void set_domain(SP_DomainDisc domain) {
            this->m_domain_disc = domain;
        }

        void set_solver(SP_Solver solver) {
            this->m_linear_solver = solver;
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif