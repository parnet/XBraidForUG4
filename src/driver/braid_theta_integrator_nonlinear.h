//todo check
#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAID_THETA_INTEGRATOR_NONLINEAR_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAID_THETA_INTEGRATOR_NONLINEAR_H

#include "lib_disc/operator/composite_conv_check.h"
#include "../integrator/theta_integrator_nl.h"
#include "lib_algebra/operator/convergence_check.h"
#include "gridfunction_base.h"

namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BraidNLIntegrator : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef ThetaIntegratorNL<TDomain, TAlgebra> T_TimeIntegrator;
        typedef SmartPtr<T_TimeIntegrator> SP_TimeIntegrator;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef StdConvCheck<typename TAlgebra::vector_type> T_ConvCheck;
        typedef SmartPtr<T_ConvCheck> SP_ConvCheck;

        //--------------------------------------------------------------------------------------------------------------

        SP_ConvCheck m_conv_check; // set but unsed here

        std::vector<SP_TimeIntegrator> time_integrator_list; //used
        SP_TimeIntegrator default_integrator = SPNULL;

        double max_error_estimate = -1;
        double m_ref_threshold = 0.2; // only set but not used

        double chk_relative = 1e-6;
        double chk_loose_tol = 1e-3;
        double chk_tight_tol = 1e-14;

        int m_ref_factor = 2; // only set but not used

        //--------------------------------------------------------------------------------------------------------------

        BraidNLIntegrator() : BraidGridFunctionBase<TDomain, TAlgebra>() {
            this->time_integrator_list = std::vector<SP_TimeIntegrator>();
        }

        BraidNLIntegrator(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidGridFunctionBase<
            TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->time_integrator_list = std::vector<SP_TimeIntegrator>();
        }

        ~BraidNLIntegrator() override= default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& pstatus) override {

            int level;
            int iteration;
            int tindex;
            int idone;
            double t_start, t_stop;
            pstatus.GetLevel(&level);
            pstatus.GetTstartTstop(&t_start, &t_stop);
            pstatus.GetDone(&idone);
            pstatus.GetTIndex(&tindex);
            pstatus.GetIter(&iteration);
            double current_dt = t_stop - t_start;
            auto* sp_u_approx_tstart = (SP_GridFunction*)u_->value;
            auto* constsp_u_approx_tstop = (SP_GridFunction*)ustop_->value;
            SP_GridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();

            SP_TimeIntegrator loc_time_integrator = this->get_integrator(level);

            if (level == 0) {
                double tol = 1.0;
                pstatus.GetSpatialAccuracy(this->chk_loose_tol, this->chk_tight_tol, &tol);
                this->set_cmp_vals(tol, chk_relative);
            }

            loc_time_integrator->set_time_step(current_dt);
            loc_time_integrator->init(*sp_u_approx_tstart->get());
            loc_time_integrator->prepare(*sp_u_approx_tstart->get());

            bool success = loc_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(),t_start);

            *sp_u_approx_tstart = sp_u_tstop_approx;
            u_->time = t_stop;

            return 0;
        };

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& pstatus) override {
            return 0;
        };

        int Sync(BraidSyncStatus& sstatus) override {
            this->m_iteration += 1;
            this->m_log->o << "Iteration: " << this->m_iteration << std::endl;

            if (this->restimate) {
                int nestimates = 0;
                sstatus.GetNumErrorEst(&nestimates);
                double* estimates = new double[nestimates];
                sstatus.GetAllErrorEst(estimates);
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
            this->m_conv_check = conv;
        }

        void set_cmp_vals(double absolute, double relative) {
            this->m_conv_check->set_minimum_defect(absolute);
            this->m_conv_check->set_reduction(relative);
            this->m_conv_check->set_maximum_steps(100);
        }


        void set_tol(double loose = 1e-3, double tight = 1e-14) {
            this->chk_loose_tol = loose;
            this->chk_tight_tol = tight;
            this->set_cmp_vals(loose, chk_relative);
        }

        void set_ref_factor(int ref_factor) {
            this->m_ref_factor = ref_factor;
        }

        void set_threshold(double ref_threshold) {
            this->m_ref_threshold = ref_threshold;
        }


        void release() {}

        void print_settings() {}

        void set_default_integrator(SP_TimeIntegrator integrator) {
            this->default_integrator = integrator;
        }

        void set_integrator(int level, SP_TimeIntegrator integrator) {
            if (this->time_integrator_list.size() < level + 1) {
                this->time_integrator_list.resize(level + 1, SPNULL);
            }
            this->time_integrator_list[level] = integrator;
        }

        SP_TimeIntegrator get_integrator(int level) {
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
#endif //UG_PLUGIN_XBRAIDFORUG4_BRAID_THETA_INTEGRATOR_NONLINEAR_H
