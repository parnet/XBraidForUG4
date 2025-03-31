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

        std::vector<SP_TimeIntegrator> time_integrator_list; //used
        SP_TimeIntegrator default_integrator = SPNULL;

        double max_error_estimate = -1;
        double m_ref_threshold = 0.2;

        double chk_relative = 1e-6;
        double chk_loose_tol = 1e-3;
        double chk_tight_tol = 1e-14;

        int coarse_linear_solver_iter = 3;
        int coarse_nonlinear_solver_iter = 3;

        //--------------------------------------------------------------------------------------------------------------

        BraidNLIntegrator() : BraidGridFunctionBase<TDomain, TAlgebra>() {
            this->time_integrator_list = std::vector<SP_TimeIntegrator>();
        }

        BraidNLIntegrator(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidGridFunctionBase<
            TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->time_integrator_list = std::vector<SP_TimeIntegrator>();
        }

        ~BraidNLIntegrator() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override {

            int level;
            status.GetLevel(&level);

            int iteration;
            status.GetIter(&iteration);

            int tindex;
            status.GetTIndex(&tindex);

            double t_start, t_stop;
            status.GetTstartTstop(&t_start, &t_stop);
            double current_dt = t_stop - t_start;

            extended_info(std::cout << "BasicDriver::Step --- "
                << " level=" << level
                << " iteration=" << iteration
                << " t_start=" << t_start
                << " t_stop=" << t_stop
                << std::endl << std::flush;)





            auto* sp_u_approx_tstart = (SP_GridFunction*)u_->value;
            auto* constsp_u_approx_tstop = (SP_GridFunction*)ustop_->value;
            SP_GridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();

            __debug(std::cout << "csp_u_tstop_approx::size = " << csp_u_tstop_approx->size() << std::endl << std::flush);
            __debug(std::cout << "sp_u_approx_tstart::size = " << sp_u_approx_tstart->size() << std::endl << std::flush);

            SP_TimeIntegrator loc_time_integrator = this->get_integrator(level);

            if (level == 0) {
                double tol = 1.0;
                status.GetSpatialAccuracy(this->chk_loose_tol, this->chk_tight_tol, &tol);
                //this->set_cmp_vals(tol, chk_relative);
                SmartPtr<IConvergenceCheck<typename TAlgebra::vector_type>> conv_check = loc_time_integrator->get_solver()->get_linear_solver()->convergence_check();
                //this->m_conv_check->set_minimum_defect(absolute);
                //this->m_conv_check->set_reduction(relative);
                //this->m_conv_check->set_maximum_steps(100);
            } else {
                SmartPtr<IConvergenceCheck<typename TAlgebra::vector_type>> conv_check = loc_time_integrator->get_solver()->get_linear_solver()->convergence_check();
                //this->m_conv_check->set_minimum_defect();
                //this->m_conv_check->set_reduction(relative);
                //this->m_conv_check->set_maximum_steps(100);
            }

            loc_time_integrator->set_time_step(current_dt);
            loc_time_integrator->init(*sp_u_approx_tstart->get());
            loc_time_integrator->prepare(*sp_u_approx_tstart->get());

            bool success = loc_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(),t_start);

            *sp_u_approx_tstart = sp_u_tstop_approx;
            //u_->time = t_stop;

            return 0;
        };

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override {
            return 0;
        };

        int Sync(BraidSyncStatus& status) override {
            this->m_iteration += 1;
            this->m_log->o << "Iteration: " << this->m_iteration << std::endl;

            if (this->restimate) {
                int nestimates = 0;
                status.GetNumErrorEst(&nestimates);
                double* estimates = new double[nestimates];
                status.GetAllErrorEst(estimates);
                double this_max = -1.0;
                for (int m = 0; m < nestimates; m++) {
                    if (this_max < estimates[m]) {
                        this_max = estimates[m];
                    }
                }
                this->m_log->o << "proc max error estimate: " << this_max << std::endl;

                MPI_Allreduce(&this_max,
                              &(this->max_error_estimate),
                              1, MPI_DOUBLE, MPI_MAX,
                              this->m_comm->TEMPORAL);
                this->m_log->o << "glob max error estimate: " << max_error_estimate << std::endl;
            }
            return 0;
        };

        //--------------------------------------------------------------------------------------------------------------

        void set_conv_check(SP_ConvCheck conv) {
            //this->m_conv_check = conv;
        }

        void set_cmp_vals(double absolute, double relative) {

        }


        void set_tol(double loose = 1e-3, double tight = 1e-14) {
            this->chk_loose_tol = loose;
            this->chk_tight_tol = tight;
            this->set_cmp_vals(loose, chk_relative);
        }



        void set_threshold(double ref_threshold) {
            this->m_ref_threshold = ref_threshold;
        }


        void release() {}

        void print_settings() {}

        void set_default_integrator(SP_TimeIntegrator integrator) {
            this->default_integrator = integrator;
        }

        void set_integrator(size_t level, SP_TimeIntegrator integrator) {
            if (this->time_integrator_list.size() < level + 1) {
                this->time_integrator_list.resize(level + 1, SPNULL);
            }
            this->time_integrator_list[level] = integrator;
        }

        SP_TimeIntegrator get_integrator(size_t level) {
            if (this->time_integrator_list.size() < level) {
                return this->default_integrator;
            } else if (this->time_integrator_list[level] != SPNULL) {
                return this->time_integrator_list[level];
            } else {
                return this->default_integrator;
            }
        }
        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif