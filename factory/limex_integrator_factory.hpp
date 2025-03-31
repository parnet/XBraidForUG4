#ifndef UGPLUGIN_XBRAIDFORUG4_FACTORY_LIMEX_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_FACTORY_LIMEX_INTEGRATOR_FACTORY_HPP

#include "Limex/time_disc/limex_integrator.hpp"
#include "interface/integrator_factory.hpp"




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class LimexFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_VectorType = typename TAlgebra::vector_type ;

        using T_TimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_ErrorEstimator = ISubDiagErrorEst<T_VectorType> ;
        using SP_ErrorEstimator = SmartPtr<T_ErrorEstimator> ;

        using T_IntegratorType = INonlinearTimeIntegrator<TDomain, TAlgebra> ;
        using T_SolverType = typename T_IntegratorType::solver_type ;
        using SP_Solver = SmartPtr<T_SolverType> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        //--------------------------------------------------------------------------------------------------------------


        SP_DomainDisc domain_disc;
        SP_Solver solver;
        SP_ErrorEstimator error_est;


        bool disable_matrix_cache = true;
        int n_stages = 4;
        double tol = 1e-3;
        double dt_min = -1;
        double dt_max = -1;
        double safety_factor = 0.25; // todo missing setter?
        double increase_factor = 2.0; // todo missing setter?
        double greedy_order_factor = 0.0; // todo missing setter?
        double level_factor = 0.5;

        std::vector<int> steps = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        //--------------------------------------------------------------------------------------------------------------

        LimexFactory() : IntegratorFactory<TDomain, TAlgebra>() {};

        ~LimexFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_domain_disc(SP_DomainDisc sp_domain_disc) {
            this->domain_disc = sp_domain_disc;
        }

        void set_solver(SP_Solver sp_solver) {
            this->solver = sp_solver;
        }

        void set_dt_min(double dt) {
            this->dt_min = dt;
        }

        void set_dt_max(double dt) {
            this->dt_max = dt;
        }

        void set_tol(double tol) {
            this->tol = tol;
        }

        void set_error_estimator(SP_ErrorEstimator sp_error_est) {
            this->error_est = sp_error_est;
        }

        void set_level_factor(double factor) {
            this->level_factor = factor;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            auto limex = make_sp(new LimexTimeIntegrator<TDomain, TAlgebra>(this->n_stages));
            for (int i = 0; i < this->n_stages; i++) {
                limex->add_stage(this->steps[i], this->solver, this->domain_disc);
            }

            limex->add_error_estimator(this->error_est);
            double factor = 1.0;
            if (level > 0) {
                factor = pow(this->level_factor, level);
            }

            double curfactor = this->tol * factor;
            limex->set_tolerance(curfactor);
            limex->set_time_step(current_dt);

            if (this->dt_min == -1) {
                limex->set_time_step(current_dt * 1e-2);
            } else {
                limex->set_time_step(this->dt_min);
            }

            if (this->dt_max == -1) {
                limex->set_time_step(current_dt);
            } else {
                limex->set_time_step(this->dt_max);
            }

            if (this->disable_matrix_cache) {
                limex->disable_matrix_cache();
            } else {
                limex->enable_matrix_cache();
            }

            limex->set_stepsize_safety_factor(this->safety_factor);
            limex->set_increase_factor(this->increase_factor);
            limex->set_stepsize_greedy_order_factor(this->greedy_order_factor);
            return limex;
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            auto limex = make_sp(new LimexTimeIntegrator<TDomain, TAlgebra>(this->n_stages));
            for (int i = 0; i < this->n_stages; i++) {
                limex->add_stage(this->steps[i], this->solver, this->domain_disc);
            }
            limex->add_error_estimator(this->error_est);
            limex->set_tolerance(this->tol);
            limex->set_time_step(current_dt);

            if (this->dt_min == -1) {
                limex->set_time_step(current_dt * 1e-2);
            } else {
                limex->set_time_step(this->dt_min);
            }

            if (this->dt_max == -1) {
                limex->set_time_step(current_dt);
            } else {
                limex->set_time_step(this->dt_max);
            }

            if (this->disable_matrix_cache) {
                limex->disable_matrix_cache();
            } else {
                limex->enable_matrix_cache();
            }

            limex->set_stepsize_safety_factor(this->safety_factor);
            limex->set_increase_factor(this->increase_factor);
            limex->set_stepsize_greedy_order_factor(this->greedy_order_factor);

            return limex;
        }
    };
}}

#endif