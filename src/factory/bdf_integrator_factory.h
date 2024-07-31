//todo check
#ifndef UG_PLUGIN_XBRAIDFORUG4_BDF_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_BDF_INTEGRATOR_FACTORY_H

#include "../interface/integrator_factory.h"
#include "../integrator/bdf_integrator.h"



namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BDF_IntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef ITimeIntegrator<TDomain, TAlgebra> T_TimeIntegrator;
        typedef SmartPtr<T_TimeIntegrator> SP_TimeIntegrator;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef LinearSolver<typename TAlgebra::vector_type> T_Solver;
        typedef SmartPtr<T_Solver> SP_Solver;

        typedef BDF_Integrator<TDomain, TAlgebra> T_BDF_Integrator;
        typedef SmartPtr<T_BDF_Integrator> SP_BDF_Integrator;

        //--------------------------------------------------------------------------------------------------------------

        SP_Solver m_linear_solver;
        SP_DomainDisc m_domain_disc;
        double order = 1;
        std::vector<double> levelorder;

        //--------------------------------------------------------------------------------------------------------------

        BDF_IntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {
            levelorder = std::vector<double>();
            for (int i = 0; i < 20; i++) { // todo remove the limit
                levelorder.template emplace_back(1);
            }
        }

        ~BDF_IntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_order(double porder) {
            this->order = porder;
        }

        void set_level_order(int level, double porder) {
            this->levelorder[level] = porder;
        }
        void set_domain(SP_DomainDisc domain) {
            this->m_domain_disc = domain;
        }
        void set_solver(SP_Solver solver) {
            this->m_linear_solver = solver;
        }
        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            SP_BDF_Integrator integrator = make_sp(new BDF_Integrator<TDomain, TAlgebra>());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_order(this->levelorder[level]);
            return integrator;
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            SP_BDF_Integrator integrator = make_sp(new BDF_Integrator<TDomain, TAlgebra>());
            integrator->set_domain(m_domain_disc);
            integrator->set_solver(m_linear_solver);
            integrator->set_order(this->order);
            return integrator;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_BDF_INTEGRATOR_FACTORY_H
