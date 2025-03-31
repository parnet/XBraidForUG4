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
            this->m_log->o << "Residual Method was called but residual support is not supported by integrator factory class" << std::flush << std::endl;
            return 1;
        };
        //--------------------------------------------------------------------------------------------------------------

        void release(){}

        void print_settings() {}

        void set_fine_time_integrator(SP_IntegratorFactory integrator) {
            this->m_fine_time_integrator_factory = integrator;
        }

        void set_coarse_time_integrator(SP_IntegratorFactory integrator) {
            this->m_coarse_time_integrator_factory = integrator;
        }// Paralog

        SP_IntegratorFactory get_fine_time_integrator() {
            return m_fine_time_integrator_factory;
        }

        SP_IntegratorFactory get_coarse_time_integrator() {
            return m_coarse_time_integrator_factory;
        }

        //--------------------------------------------------------------------------------------------------------------

        SP_IntegratorFactory m_fine_time_integrator_factory;
        SP_IntegratorFactory m_coarse_time_integrator_factory;

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

        SP_TimeIntegrator loc_time_integrator;
        if (level <= 0) {
            loc_time_integrator = m_fine_time_integrator_factory->create_level_time_integrator( current_dt, bool(idone), level);
        } else {
            loc_time_integrator = m_coarse_time_integrator_factory->create_level_time_integrator( current_dt, bool(idone), level);
        }
        auto* sp_u_approx_tstart = (SP_GridFunction*)u_->value;
        auto* constsp_u_approx_tstop = (SP_GridFunction*)ustop_->value;
        SP_GridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        if (fstop_ != nullptr) {
            this->m_log->o << "Warning residual is ignored" << std::endl << std::flush;
            exit(127);
        }
        loc_time_integrator->init(*sp_u_approx_tstart->get()->clone());
        loc_time_integrator->prepare(*sp_u_approx_tstart->get()->clone());
        bool success = loc_time_integrator->apply(sp_u_tstop_approx, t_stop, sp_u_approx_tstart->cast_const(),t_start);
        if (!success) {
            this->m_log->o << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        }
        *sp_u_approx_tstart = sp_u_tstop_approx;
        return 0;
    }
}}
#endif