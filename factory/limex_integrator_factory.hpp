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


        LimexFactory() : IntegratorFactory<TDomain, TAlgebra>() {};

        ~LimexFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_domain(SP_DomainDisc sp_domain_disc) {
            this->domain_disc_ = sp_domain_disc;
        }

        void set_solver(SP_Solver sp_solver) {
            this->solver_ = sp_solver;
        }

        void set_dt_min(double dt) {
            this->dt_min_ = dt;
        }

        void set_dt_max(double dt) {
            this->dt_max_ = dt;
        }

        void set_tol(double tol) {
            this->tol_ = tol;
        }

        void set_error_estimator(SP_ErrorEstimator sp_error_est) {
            this->error_est_ = sp_error_est;
        }

        void set_level_factor(double factor) {
            this->level_factor_ = factor;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            auto limex = make_sp(new LimexTimeIntegrator<TDomain, TAlgebra>(this->n_stages_));
            for (int i = 0; i < this->n_stages_; i++) {
                limex->add_stage(this->steps_[i], this->solver_, this->domain_disc_);
            }

            limex->add_error_estimator(this->error_est_);
            double factor = 1.0;
            if (level > 0) {
                factor = pow(this->level_factor_, level);
            }

            double curfactor = this->tol_ * factor;
            limex->set_tolerance(curfactor);
            limex->set_time_step(current_dt);

            if (this->dt_min_ == -1) {
                limex->set_time_step(current_dt * 1e-2);
            } else {
                limex->set_time_step(this->dt_min_);
            }

            if (this->dt_max_ == -1) {
                limex->set_time_step(current_dt);
            } else {
                limex->set_time_step(this->dt_max_);
            }

            if (this->disable_matrix_cache_) {
                limex->disable_matrix_cache();
            } else {
                limex->enable_matrix_cache();
            }

            limex->set_stepsize_safety_factor(this->safety_factor_);
            limex->set_increase_factor(this->increase_factor_);
            limex->set_stepsize_greedy_order_factor(this->greedy_order_factor_);
            return limex;
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            auto limex = make_sp(new LimexTimeIntegrator<TDomain, TAlgebra>(this->n_stages_));
            for (int i = 0; i < this->n_stages_; i++) {
                limex->add_stage(this->steps_[i], this->solver_, this->domain_disc_);
            }
            limex->add_error_estimator(this->error_est_);
            limex->set_tolerance(this->tol_);
            limex->set_time_step(current_dt);

            if (this->dt_min_ == -1) {
                limex->set_time_step(current_dt * 1e-2);
            } else {
                limex->set_time_step(this->dt_min_);
            }

            if (this->dt_max_ == -1) {
                limex->set_time_step(current_dt);
            } else {
                limex->set_time_step(this->dt_max_);
            }

            if (this->disable_matrix_cache_) {
                limex->disable_matrix_cache();
            } else {
                limex->enable_matrix_cache();
            }

            limex->set_stepsize_safety_factor(this->safety_factor_);
            limex->set_increase_factor(this->increase_factor_);
            limex->set_stepsize_greedy_order_factor(this->greedy_order_factor_);

            return limex;
        }
        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc domain_disc_;
        SP_Solver solver_;
        SP_ErrorEstimator error_est_;

        bool disable_matrix_cache_ = true;
        int n_stages_ = 4;
        double tol_ = 1e-3;
        double dt_min_ = -1;
        double dt_max_ = -1;
        double safety_factor_ = 0.25;
        double increase_factor_ = 2.0;
        double greedy_order_factor_ = 0.0;
        double level_factor_ = 0.5;

        std::vector<int> steps_ = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif