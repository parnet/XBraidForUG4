#ifndef UGPLUGIN_XBRAIDFORUG4_FACTORY_BDF_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_FACTORY_BDF_INTEGRATOR_FACTORY_HPP

#include "interface/integrator_factory.hpp"
#include "integrator/bdf_integrator.hpp"




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BDF_IntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_TimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_BDF_Integrator = BDF_Integrator<TDomain, TAlgebra> ;
        using SP_BDF_Integrator = SmartPtr<T_BDF_Integrator> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_Solver m_linear_solver;
        SP_DomainDisc m_domain_disc;
        double order = 1;
        std::vector<double> levelorder;

        //--------------------------------------------------------------------------------------------------------------

        BDF_IntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {
            levelorder = std::vector<double>();
            for (int i = 0; i < 20; i++) { // todo remove the limit
                levelorder.emplace_back(1);
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
#endif