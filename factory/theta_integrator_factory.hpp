#ifndef UGPLUGIN_XBRAIDFORUG4_FACTORY_THETA_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_FACTORY_THETA_INTEGRATOR_FACTORY_HPP

#include "interface/integrator_factory.hpp"
#include "integrator/theta_single_timestep.hpp"




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class ThetaIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_TimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_ThetaIntegrator = ThetaSingleTimeStep<TDomain, TAlgebra> ;
        using SP_ThetaIntegrator = SmartPtr<T_ThetaIntegrator> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_Solver m_linear_solver;
        SP_DomainDisc m_domain_disc;

        int max_level = 20;

        double default_theta = 1;
        std::vector<double> leveltheta;

        //--------------------------------------------------------------------------------------------------------------

        ThetaIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {
            leveltheta = std::vector<double>();
            for (int i = 0; i < max_level; i++) { // todo remove the limit
                leveltheta.emplace_back(1);
            }
        }
        ~ThetaIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_theta(double theta) {
            this->default_theta = theta;
        }

        void set_max_level(int max_level) {
            this->max_level = max_level;
            // todo expand structures if >
        }

        void set_level_theta(int level, double theta) {
            this->leveltheta[level] = theta;
        }

        void set_domain(SP_DomainDisc domain) {
            this->m_domain_disc = domain;
        }

        void set_solver(SP_Solver solver) {
            this->m_linear_solver = solver;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_theta(this->leveltheta[level]);
            return integrator;
        }


        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_theta(this->default_theta);
            return integrator;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif