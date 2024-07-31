#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAID_INTEGRATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAID_INTEGRATOR_H


#include "gridfunction_base.h"


namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BraidIntegrator : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef typename TAlgebra::matrix_type T_Matrix;
        typedef typename TAlgebra::vector_type T_Vector;
        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;
        typedef ITimeIntegrator<TDomain, TAlgebra> T_TimeIntegrator;
        typedef SmartPtr<T_TimeIntegrator> SP_TimeIntegrator;
        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SPDomainDisc;
        typedef StdConvCheck<typename TAlgebra::vector_type> T_ConvCheck;
        typedef SmartPtr<T_ConvCheck> SP_ConvCheck;

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_TimeIntegrator> time_integrator_list;
        SP_TimeIntegrator default_integrator;

        double max_error_estimate = -1; // error estimation and refine
        double m_ref_threshold = 0.2;// error estimation and refine
        int m_ref_factor = 2;// error estimation and refine

        //--------------------------------------------------------------------------------------------------------------

        BraidIntegrator() : BraidGridFunctionBase<TDomain, TAlgebra>() {
            this->time_integrator_list = std::vector<SP_TimeIntegrator>();
        }
        BraidIntegrator(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase< TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {

            this->time_integrator_list = std::vector<SP_TimeIntegrator>();
        }
        ~BraidIntegrator() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& pstatus) override {
            int level; // level
            pstatus.GetLevel(&level);
            double t_start, t_stop;
            pstatus.GetTstartTstop(&t_start, &t_stop);
            double current_dt = t_stop - t_start;
            int idone;
            pstatus.GetDone(&idone);
            int tindex;
            pstatus.GetTIndex(&tindex);
            int iteration;
            pstatus.GetIter(&iteration);
            auto* sp_u_approx_tstart = (SP_GridFunction*)u_->value;
            auto* constsp_u_approx_tstop = (SP_GridFunction*)ustop_->value;
            SP_GridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
            SP_TimeIntegrator loc_time_integrator = this->get_integrator(level);
            loc_time_integrator->set_time_step(current_dt);
            loc_time_integrator->init(*sp_u_approx_tstart->get());
            loc_time_integrator->prepare(*sp_u_approx_tstart->get());
            bool success = loc_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(), t_start);
            if (!success) {
                //&& m_force_conv) {
                this->m_log->o << "!!! Failure convergence not reached" << std::endl;
                exit(127);
            }
            *sp_u_approx_tstart = sp_u_tstop_approx;
            u_->time = t_stop;
            if (this->restimate) {
                double local_estimate;
                pstatus.GetSingleErrorEstStep(&local_estimate);
                if (local_estimate != -1.0) {
                    this->m_log->o << ">> estim idx=" << tindex + 1 << " t=" << t_stop << " iter=" << iteration <<
                        " level=" << level << " done=" << idone << " value=" << local_estimate << std::endl;
                }
                if (this->refine_time && level == 0 && iteration > 1 && local_estimate > m_ref_threshold * this->
                    max_error_estimate && local_estimate != -1.0) {
                    this->m_log->o << "refine by " << m_ref_factor << std::endl;
                    pstatus.SetRFactor(m_ref_factor);
                }
            }
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
            }

            if (this->time_integrator_list[level] != SPNULL) {
                return this->time_integrator_list[level];
            }

            return this->default_integrator;
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_BRAID_INTEGRATOR_H
