#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_INTEGRATOR_FACTORY_HPP


#include "Limex/time_disc/time_integrator.hpp"

#include "interface/integrator_factory.hpp"
#include "core/space_time_communicator.hpp"

#include "gridfunction_base.hpp"




namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidIntegratorFactory : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using SP_Communicator = SmartPtr<SpaceTimeCommunicator> ;

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_IntegratorFactory = IntegratorFactory<TDomain, TAlgebra> ;
        using SP_IntegratorFactory = SmartPtr<T_IntegratorFactory> ;

        using T_TimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;
        using LevelIntegratorList = std::vector<SP_TimeIntegrator> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;


        //--------------------------------------------------------------------------------------------------------------

        BraidIntegratorFactory() : BraidGridFunctionBase<TDomain, TAlgebra>() {}

        BraidIntegratorFactory(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {}

        ~BraidIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override;;

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override {
            this->log_->o << "Residual Method was called but residual support is not supported by integrator factory class" << std::flush << std::endl;
            return 1;
        };
        //--------------------------------------------------------------------------------------------------------------

        void release(){}

        void print_settings() {}

        void set_default_integrator(SP_IntegratorFactory integrator) {
            this->default_integrator_factory_ = integrator;
        }

        void set_integrator(int level, SP_IntegratorFactory integrator) {
            if (this->integrator_factory_list_.size() < level + 1) {
                this->integrator_factory_list_.resize(level + 1, SPNULL);
            }
            this->integrator_factory_list_[level] = integrator;
        }

        //2025-04 void set_fine_time_integrator(SP_IntegratorFactory integrator) {
        //2025-04     this->fine_time_integrator_factory_ = integrator;
        //2025-04 }

        //2025-04 void set_coarse_time_integrator(SP_IntegratorFactory integrator) {
        //2025-04     this->coarse_time_integrator_factory_ = integrator;
        //2025-04 }

        //2025-04 SP_IntegratorFactory get_fine_time_integrator() {
        //2025-04     return fine_time_integrator_factory_;
        //2025-04 }

        //2025-04 SP_IntegratorFactory get_coarse_time_integrator() {
        //2025-04     return coarse_time_integrator_factory_;
        //2025-04 }

        SP_IntegratorFactory get_integrator(int level) {
            if (this->integrator_factory_list_.size() < level) {
                return this->default_integrator_factory_;
            }

            if (this->integrator_factory_list_[level] != SPNULL) {
                return this->integrator_factory_list_[level];
            }

            return this->default_integrator_factory_;
        }

        SP_IntegratorFactory get_default_integrator() {
            return default_integrator_factory_;
        }

        //--------------------------------------------------------------------------------------------------------------

        //2025-04 SP_IntegratorFactory fine_time_integrator_factory_;
        //2025-04 SP_IntegratorFactory coarse_time_integrator_factory_;

        SP_IntegratorFactory default_integrator_factory_;
        std::vector<SP_IntegratorFactory> integrator_factory_list_;

        //--------------------------------------------------------------------------------------------------------------
    };

    template<typename TDomain, typename TAlgebra>
    int BraidIntegratorFactory<TDomain, TAlgebra>::Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_,
        BraidStepStatus &status) {
        int level;
        status.GetLevel(&level);

        double t_start, t_stop;
        status.GetTstartTstop(&t_start, &t_stop);
        double current_dt = t_stop - t_start;

        int idone;
        status.GetDone(&idone);

        int tindex;
        status.GetTIndex(&tindex);

        int iteration;
        status.GetIter(&iteration);

        SP_TimeIntegrator loc_time_integrator = this->get_integrator(level)->create_level_time_integrator(current_dt, bool(idone), level);
        //2025-04 if (level <= 0) {
        //2025-04     loc_time_integrator = fine_time_integrator_factory_->create_level_time_integrator( current_dt, bool(idone), level);
        //2025-04 } else {
        //2025-04     loc_time_integrator = coarse_time_integrator_factory_->create_level_time_integrator( current_dt, bool(idone), level);
        //2025-04 }
        auto* sp_u_approx_tstart = (SP_GridFunction*)u_->value_;
        auto* constsp_u_approx_tstop = (SP_GridFunction*)ustop_->value_;
        SP_GridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        if (fstop_ != nullptr) {
            this->log_->o << "Warning residual is ignored" << std::endl << std::flush;
            exit(127);
        }
        loc_time_integrator->init(*sp_u_approx_tstart->get()->clone());
        loc_time_integrator->prepare(*sp_u_approx_tstart->get()->clone());
        bool success = loc_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(),t_start);
        if (!success) {
            this->log_->o << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        }
        *sp_u_approx_tstart = sp_u_tstop_approx;
        return 0;
    }
}}
#endif