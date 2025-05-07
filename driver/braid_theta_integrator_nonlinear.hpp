#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_THETA_INTEGRATOR_NONLINEAR_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_THETA_INTEGRATOR_NONLINEAR_HPP

#include "lib_disc/operator/composite_conv_check.h"
#include "integrator/theta_integrator_nl.hpp"
#include "lib_algebra/operator/convergence_check.h"
#include "gridfunction_base.hpp"


namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidNLIntegrator : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_TimeIntegrator = ThetaIntegratorNL<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_ConvCheck = StdConvCheck<typename TAlgebra::vector_type> ;
        using SP_ConvCheck = SmartPtr<T_ConvCheck> ;


        //--------------------------------------------------------------------------------------------------------------

        BraidNLIntegrator() : BraidGridFunctionBase<TDomain, TAlgebra>() {
            this->time_integrator_list_ = std::vector<SP_TimeIntegrator>();
        }

        BraidNLIntegrator(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidGridFunctionBase<
            TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->time_integrator_list_ = std::vector<SP_TimeIntegrator>();
        }

        ~BraidNLIntegrator() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override;

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override {
            return 0;
        };

        int Sync(BraidSyncStatus& status) override;;

        //--------------------------------------------------------------------------------------------------------------

        void set_conv_check(SP_ConvCheck conv) {
        }

        void set_cmp_vals(double absolute, double relative) {

        }


        void set_tol(double loose = 1e-3, double tight = 1e-14) {
            this->chk_loose_tol_ = loose;
            this->chk_tight_tol_ = tight;
            this->set_cmp_vals(loose, chk_relative_);
        }



        void set_threshold(double ref_threshold) {
            this->ref_threshold_ = ref_threshold;
        }


        void release() {}

        void print_settings() {}

        void set_default_integrator(SP_TimeIntegrator integrator) {
            this->default_integrator_ = integrator;
        }

        void set_integrator(size_t level, SP_TimeIntegrator integrator);

        SP_TimeIntegrator get_integrator(size_t level);

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_TimeIntegrator> time_integrator_list_; //used
        SP_TimeIntegrator default_integrator_ = SPNULL;

        double max_error_estimate_ = -1; //æ concept [refinement]
        double ref_threshold_ = 0.2; //æ concept [refinement]

        double chk_relative_ = 1e-6; //æ concept [adaptive_tolerance]
        double chk_loose_tol_ = 1e-3; //æ concept [adaptive_tolerance]
        double chk_tight_tol_ = 1e-14;//æ concept [adaptive_tolerance]

        int coarse_linear_solver_iter_ = 3; //æ concept [adaptive_tolerance]
        int coarse_nonlinear_solver_iter_ = 3;//æ concept [adaptive_tolerance]

    };

    template<typename TDomain, typename TAlgebra>
    int BraidNLIntegrator<TDomain, TAlgebra>::Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
        BraidStepStatus &status) {

        int level;
        status.GetLevel(&level);

        int iteration;
        status.GetIter(&iteration);

        int t_index;
        status.GetTIndex(&t_index);

        double t_start, t_stop;
        status.GetTstartTstop(&t_start, &t_stop);
        double current_dt = t_stop - t_start;

        extended_info(std::cout << "BasicDriver::Step --- "
            << " level=" << level
            << " iteration=" << iteration
            << " t_start=" << t_start
            << " t_stop=" << t_stop
            << std::endl << std::flush;)


        //2025-05
        //std::cout << "<debug_section>" << std::endl;
        //double u_start = 0.0;
        //double u_stop = 0.0;
        //this->SpatialNorm(u_, &u_start);
        //this->SpatialNorm(ustop_, &u_stop);
        //std::cout << "norm(u_start)=" << u_start << std::endl;
        //std::cout << "norm(u_stop)=" << u_stop << std::endl;
        //std::cout << "t_start=" << t_start << std::endl;
        //std::cout << "t_stop=" << t_stop << std::endl;
        //std::cout << "dt=" << current_dt << std::endl;
        //std::cout << "</debug_section>" << std::endl;



        auto* sp_u_approx_tstart = (SP_GridFunction*)u_->value_;
        auto* constsp_u_approx_tstop = (SP_GridFunction*)ustop_->value_;
        SP_GridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();

        __debug(std::cout << "csp_u_tstop_approx::size = " << csp_u_tstop_approx->size() << std::endl << std::flush);
        __debug(std::cout << "sp_u_approx_tstart::size = " << sp_u_approx_tstart->size() << std::endl << std::flush);

        SP_TimeIntegrator loc_time_integrator = this->get_integrator(level);

        //2025-05
        //if (level == 0) {
        //    double tol = 1.0;
        //    status.GetSpatialAccuracy(this->chk_loose_tol_, this->chk_tight_tol_, &tol);
            //this->set_cmp_vals(tol, chk_relative);
        //    SmartPtr<IConvergenceCheck<typename TAlgebra::vector_type>> conv_check = loc_time_integrator->get_solver()->get_linear_solver()->convergence_check();
            //this->m_conv_check->set_minimum_defect(absolute);
            //this->m_conv_check->set_reduction(relative);
            //this->m_conv_check->set_maximum_steps(100);
        //} else {
        //    SmartPtr<IConvergenceCheck<typename TAlgebra::vector_type>> conv_check = loc_time_integrator->get_solver()->get_linear_solver()->convergence_check();
            //this->m_conv_check->set_minimum_defect();
            //this->m_conv_check->set_reduction(relative);
            //this->m_conv_check->set_maximum_steps(100);
        //}

        loc_time_integrator->set_time_step(current_dt);
        loc_time_integrator->init(*sp_u_approx_tstart->get()); // todo method is empty
        loc_time_integrator->prepare(*sp_u_approx_tstart->get()); // todo method is empty

        loc_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(),t_start);

        *sp_u_approx_tstart = sp_u_tstop_approx;
        //u_->time = t_stop;

        //2025-05
        //std::cout << "<debug_section>" << std::endl;
        //this->SpatialNorm(u_, &u_start);
        //this->SpatialNorm(ustop_, &u_stop);
        //std::cout << "norm(u_start)=" << u_start << std::endl;
        //std::cout << "norm(u_stop)=" << u_stop << std::endl;
        //std::cout << "</debug_section>" << std::endl;
        write_script(this->script_->Step(u_, ustop_, fstop_,status);)
        return 0;
    }

    template<typename TDomain, typename TAlgebra>
    int BraidNLIntegrator<TDomain, TAlgebra>::Sync(BraidSyncStatus &status) {
        this->iteration_ += 1;
        this->log_->o << "Iteration: " << this->iteration_ << std::endl;

        if (this->restimate_) {
            int nestimates = 0;
            status.GetNumErrorEst(&nestimates);
            auto* estimates = new double[nestimates];
            status.GetAllErrorEst(estimates);
            double this_max = -1.0;
            for (int m = 0; m < nestimates; m++) {
                if (this_max < estimates[m]) {
                    this_max = estimates[m];
                }
            }
            this->log_->o << "proc max error estimate: " << this_max << std::endl;

            MPI_Allreduce(&this_max,
                          &(this->max_error_estimate_),
                          1, MPI_DOUBLE, MPI_MAX,
                          this->comm_->TEMPORAL);
            this->log_->o << "glob max error estimate: " << max_error_estimate_ << std::endl;
        }
        return 0;
    }

    template<typename TDomain, typename TAlgebra>
    void BraidNLIntegrator<TDomain, TAlgebra>::set_integrator(size_t level, SP_TimeIntegrator integrator) {
        if (this->time_integrator_list_.size() < level + 1) {
            this->time_integrator_list_.resize(level + 1, SPNULL);
        }
        this->time_integrator_list_[level] = integrator;
    }

    template<typename TDomain, typename TAlgebra>
    typename BraidNLIntegrator<TDomain, TAlgebra>::SP_TimeIntegrator BraidNLIntegrator<TDomain, TAlgebra>::
    get_integrator(size_t level) {
        std::cout << "lst_size = " << this->time_integrator_list_.size() << std::endl << std::flush;
        std::cout << "level="<< level  << std::endl << std::flush;
        if (this->time_integrator_list_.size() <= level) {
            return this->default_integrator_;
        } else if (this->time_integrator_list_[level] != SPNULL) {
            return this->time_integrator_list_[level];
        } else {
            return this->default_integrator_;
        }
    }
}}
#endif