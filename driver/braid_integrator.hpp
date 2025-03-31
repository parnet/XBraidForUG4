#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_INTEGRATOR_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_INTEGRATOR_HPP


#include "gridfunction_base.hpp"




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidIntegrator : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_Matrix = typename TAlgebra::matrix_type ;
        using T_Vector = typename TAlgebra::vector_type ;

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_TimeIntegrator =  ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_DomainDisc= IDomainDiscretization<TAlgebra> ;
        using SPDomainDisc = SmartPtr<T_DomainDisc> ;

        using T_ConvCheck = StdConvCheck<typename TAlgebra::vector_type> ;
        using SP_ConvCheck = SmartPtr<T_ConvCheck> ;

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

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override;;

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override {
            return 0;
        };

        int Sync(BraidSyncStatus& status) override;;

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

        void set_integrator(size_t level, SP_TimeIntegrator integrator) {
            if (this->time_integrator_list.size() < level + 1) {
                this->time_integrator_list.resize(level + 1, SPNULL);
            }
            this->time_integrator_list[level] = integrator;
        }

        SP_TimeIntegrator get_integrator(size_t level) {
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

    template<typename TDomain, typename TAlgebra>
    int BraidIntegrator<TDomain, TAlgebra>::Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
        BraidStepStatus &status) {
        int level;
        status.GetLevel(&level);

        double t_start, t_stop;
        status.GetTstartTstop(&t_start, &t_stop);
        double current_dt = t_stop - t_start;


        int tindex;
        status.GetTIndex(&tindex);

        int iteration;
        status.GetIter(&iteration);

        auto* sp_u_approx_tstart = static_cast<SP_GridFunction *>(u_->value);
        auto* constsp_u_approx_tstop = static_cast<SP_GridFunction *>(ustop_->value);
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

        if (this->restimate) {
            double local_estimate;
            status.GetSingleErrorEstStep(&local_estimate);
            if (local_estimate != -1.0) {
                this->m_log->o << ">> estim idx=" << tindex + 1 << " t=" << t_stop << " iter=" << iteration <<
                        " level=" << level << " value=" << local_estimate << std::endl;
            }
            if (this->refine_time && level == 0 && iteration > 1 && local_estimate > m_ref_threshold * this->
                max_error_estimate && local_estimate != -1.0) {
                this->m_log->o << "refine by " << m_ref_factor << std::endl;
                status.SetRFactor(m_ref_factor);
            }
        }
        return 0;
    }

    template<typename TDomain, typename TAlgebra>
    int BraidIntegrator<TDomain, TAlgebra>::Sync(BraidSyncStatus &status) {
        this->m_iteration += 1;
        this->m_log->o << "Iteration: " << this->m_iteration << std::endl;
        if (this->restimate) {
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
            this->m_log->o << "proc max error estimate: " << this_max << std::endl;
            MPI_Allreduce(&this_max,
                          &(this->max_error_estimate),
                          1, MPI_DOUBLE, MPI_MAX,
                          this->m_comm->TEMPORAL);
            this->m_log->o << "glob max error estimate: " << max_error_estimate << std::endl;
            delete[] estimates;
        }
        return 0;
    }
}}
#endif