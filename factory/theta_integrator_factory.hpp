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

        ThetaIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {
            leveltheta_ = std::vector<double>();
            for (int i = 0; i < max_level_; i++) {
                leveltheta_.emplace_back(1);
            }
        }
        ~ThetaIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_theta(double theta) {
            this->default_theta_ = theta;
        }

        void set_max_level(int max_level) {
            this->max_level_ = max_level;
        }

        void set_level_theta(int level, double theta) {
            this->leveltheta_[level] = theta;
        }

        void set_domain(SP_DomainDisc domain) {
            this->domain_disc_ = domain;
        }

        void set_solver(SP_Solver solver) {
            this->linear_solver_ = solver;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(domain_disc_);
            integrator->set_solver(linear_solver_);
            integrator->set_theta(this->leveltheta_[level]);
            return integrator;
        }


        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            SP_ThetaIntegrator integrator = make_sp(new T_ThetaIntegrator());
            integrator->set_domain(domain_disc_);
            integrator->set_solver(linear_solver_);
            integrator->set_theta(this->default_theta_);
            return integrator;
        };
        //--------------------------------------------------------------------------------------------------------------

        SP_Solver linear_solver_;
        SP_DomainDisc domain_disc_;

        int max_level_ = 20;

        double default_theta_ = 1;
        std::vector<double> leveltheta_;

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif