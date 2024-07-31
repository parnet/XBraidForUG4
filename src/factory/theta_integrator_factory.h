#ifndef UG_PLUGIN_XBRAIDFORUG4_THETA_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_THETA_INTEGRATOR_FACTORY_H

#include "../interface/integrator_factory.h"
#include "../integrator/theta_integrator.h"


namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class ThetaIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef ITimeIntegrator<TDomain, TAlgebra> T_TimeIntegrator;
        typedef SmartPtr<T_TimeIntegrator> SP_TimeIntegrator;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef LinearSolver<typename TAlgebra::vector_type> T_Solver;
        typedef SmartPtr<T_Solver> SP_Solver;

        typedef ThetaIntegrator<TDomain, TAlgebra> T_ThetaIntegrator;
        typedef SmartPtr<T_ThetaIntegrator> SP_ThetaIntegrator;

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
                leveltheta.template emplace_back(1);
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
            SP_ThetaIntegrator integrator = make_sp(new ThetaIntegrator<TDomain, TAlgebra>());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_theta(this->leveltheta[level]);
            return integrator;
        }


        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            SP_ThetaIntegrator integrator = make_sp(new ThetaIntegrator<TDomain, TAlgebra>());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_theta(this->default_theta);
            return integrator;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif //UG_PLUGIN_XBRAIDFORUG4_THETA_INTEGRATOR_FACTORY_H
