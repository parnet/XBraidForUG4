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

        BDF_IntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {
            levelorder_ = std::vector<double>();
            for (int i = 0; i < 20; i++) {
                levelorder_.emplace_back(1);
            }
        }

        ~BDF_IntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_order(double order) {
            this->order_ = order;
        }

        void set_level_order(int level, double porder) {
            this->levelorder_[level] = porder;
        }
        void set_domain(SP_DomainDisc domain) {
            this->domain_disc_ = domain;
        }
        void set_solver(SP_Solver solver) {
            this->linear_solver_ = solver;
        }
        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            SP_BDF_Integrator integrator = make_sp(new BDF_Integrator<TDomain, TAlgebra>());
            integrator->set_domain(domain_disc_);
            integrator->set_solver(linear_solver_);
            integrator->set_order(this->levelorder_[level]);
            return integrator;
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            SP_BDF_Integrator integrator = make_sp(new BDF_Integrator<TDomain, TAlgebra>());
            integrator->set_domain(domain_disc_);
            integrator->set_solver(linear_solver_);
            integrator->set_order(this->order_);
            return integrator;
        };

        //--------------------------------------------------------------------------------------------------------------

        SP_Solver linear_solver_;
        SP_DomainDisc domain_disc_;
        double order_ = 1;
        std::vector<double> levelorder_;
        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif