#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAID_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAID_INTEGRATOR_FACTORY_H


#include "../../../Limex/time_disc/time_integrator.hpp"

#include "../interface/integrator_factory.h"
#include "../core/space_time_communicator.h"

#include "gridfunction_base.h"


namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BraidIntegratorFactory : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef SmartPtr<SpaceTimeCommunicator> SP_Communicator;

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef IntegratorFactory<TDomain, TAlgebra> T_IntegratorFactory;
        typedef SmartPtr<T_IntegratorFactory> SP_IntegratorFactory;

        typedef ITimeIntegrator<TDomain, TAlgebra> T_TimeIntegrator;
        typedef SmartPtr<T_TimeIntegrator> SP_TimeIntegrator;
        typedef std::vector<SP_TimeIntegrator> LevelIntegratorList;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        //--------------------------------------------------------------------------------------------------------------

        SP_IntegratorFactory m_fine_time_integrator_factory;
        SP_IntegratorFactory m_coarse_time_integrator_factory;

        //--------------------------------------------------------------------------------------------------------------

        BraidIntegratorFactory() : BraidGridFunctionBase<TDomain, TAlgebra>() {}

        BraidIntegratorFactory(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {}

        ~BraidIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& pstatus) override {
            int level;
            pstatus.GetLevel(&level);
            double t_start, t_stop;
            pstatus.GetTstartTstop(&t_start, &t_stop);
            int idone;
            pstatus.GetDone(&idone);
            double current_dt = t_stop - t_start;
            int tindex;
            pstatus.GetTIndex(&tindex);
            int iteration;
            pstatus.GetIter(&iteration);
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
            u_->time = t_stop;
            return 0;
        };

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& pstatus) override {
            this->m_log->o << "Residual Method was called but residual support is not supported by integrator factory class" << std::flush << std::endl;
            return 1;
        };
        //--------------------------------------------------------------------------------------------------------------

        void release(){}

        void print_settings() {}

        void set_max_levels(int level) { }

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
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_BRAID_INTEGRATOR_FACTORY_H
