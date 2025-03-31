//todo check || composite conv check => remove
#ifndef UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_RESIDUAL_STEPPER_HPP
#define UGPLUGIN_XBRAIDFORUG4_DRIVER_BRAID_RESIDUAL_STEPPER_HPP


//2025-03 #define include_plugin "Limex" file
//2025-03 #define include_plugin  ../../"Limex"/file
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/operator/composite_conv_check.h"
#include "gridfunction_base.hpp"



//2025-03 #include "../../Limex/time_disc/interface/nonlinear_time_integrator.h"
//2025-03 #include "../../Limex/time_disc/timestep/linear_implicit_timestep.h"
#include "Limex/time_disc/linear_implicit_timestep.h"
//2025-03 #include "../../Limex/time_disc/nonlinear_integrator/simple_integrator.hpp"
#include "Limex/time_disc/simple_integrator.hpp"
//2025-03 #include "../../Limex/time_disc/extrapolation/aitken_neville_timex.h"


#include "core/space_time_communicator.hpp"
#include "util/parallel_logger.hpp"

/* 2025-03-25
namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidResidualStepper : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_Vector =typename TAlgebra::vector_type ;

        using T_ErrorEstimatorType = ISubDiagErrorEst<T_Vector> ;
        using SP_ErrorEstimatorType = SmartPtr<T_ErrorEstimatorType> ;

        using T_TimexType = AitkenNevilleTimex<T_Vector> ;
        using SP_TimexType = SmartPtr<T_TimexType> ;

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_TimeStep = ThetaTimeStep<TAlgebra> ;
        using SP_TimeStep = SmartPtr<T_TimeStep> ;

        using T_INonlinearTimeIntegrator = INonlinearTimeIntegrator<TDomain, TAlgebra> ;
        using T_SolverType = typename T_INonlinearTimeIntegrator::solver_type ;
        using SP_SolverType = SmartPtr<T_SolverType> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_TimeIntegratorType = SimpleTimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegratorType = SmartPtr<T_TimeIntegratorType> ;

        using T_TimeStepType = LinearImplicitEuler<TAlgebra> ;
        using SP_TimeStepType = SmartPtr<T_TimeStepType> ;

        using T_Conv =  StdConvCheck<typename TAlgebra::vector_type> ;
        using SP_Conv = SmartPtr<T_Conv> ;

        using SP_SpaceTimeCommunicator = SmartPtr<SpaceTimeCommunicator> ;

        using SP_Paralog = SmartPtr<ParallelLogger> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_ApproxSpace = ApproximationSpace<TDomain> ;
        using SP_ApproxSpace = SmartPtr<T_ApproxSpace> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_TimeStep m_default_time_step;
        SP_Solver linSolver;
        SP_ApproxSpace m_approx_space; // only used for conv check, could be deleted

        //--------------------------------------------------------------------------------------------------------------

        BraidResidualStepper() : BraidGridFunctionBase<TDomain, TAlgebra>() {
            this->can_residual_method = true;
        }

        BraidResidualStepper(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
            : BraidGridFunctionBase<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {
            this->provide_residual = true;
        }

        ~BraidResidualStepper() override = default;

        //--------------------------------------------------------------------------------------------------------------


        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& status) override {
            //print_status(this->m_log->o,pstatus);
            int level;
            status.GetLevel(&level);

            int idone;
            status.GetDone(&idone);

            int tindex;
            status.GetTIndex(&tindex);

            int iteration;
            status.GetIter(&iteration);

            double t_start,  t_stop;
            status.GetTstartTstop(&t_start, &t_stop);
            double current_dt = t_stop - t_start;


            auto* sp_u_approx_tstart = reinterpret_cast<SP_GridFunction *>(u_->value);
            auto* sp_u_tstop_approx = reinterpret_cast<SP_GridFunction *>(ustop_->value);
            SP_GridFunction csp_u_tstop_approx = sp_u_tstop_approx->get()->clone();
            SP_GridFunction sp_rhs = this->m_u0->clone_without_values();

            bool success;
            this->m_default_time_step = make_sp(new ThetaTimeStep<TAlgebra>(this->m_domain_disc));
            this->m_default_time_step->set_theta(1.0); // implicit euler;

            auto solTimeSeries = make_sp(new VectorTimeSeries<typename TAlgebra::vector_type>());
            solTimeSeries->push(sp_u_approx_tstart->get()->clone(), t_start);
            this->m_default_time_step->prepare_step(solTimeSeries, current_dt);

            const GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();
            auto Operator_A = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_default_time_step, gridlevel));
            auto* ptr_Operator_A = Operator_A.get();

            this->m_default_time_step->assemble_linear(*ptr_Operator_A, *sp_rhs.get(), gridlevel);

            int dim = 2;
            double p0ablute = 1.0;
            double p0relative = 1.0;
            p0ablute = 1e-14;
            p0relative = 1e-6;

            auto cmpConvCheckC = make_sp(new CompositeConvCheck<typename TAlgebra::vector_type, TDomain>(this->m_approx_space)); // todo remove and replace by general approch
            cmpConvCheckC->set_component_check("ux", p0ablute, p0relative);
            cmpConvCheckC->set_component_check("uy", p0ablute, p0relative);
            if (dim == 3) { cmpConvCheckC->set_component_check("uz", p0ablute, p0relative); }
            cmpConvCheckC->set_component_check("p", p0ablute, p0relative);
            cmpConvCheckC->set_maximum_steps(100);
            cmpConvCheckC->set_supress_unsuccessful(true);
            linSolver->set_convergence_check(cmpConvCheckC);

            linSolver->init(Operator_A, *sp_u_approx_tstart->get());
            success = linSolver->apply(*csp_u_tstop_approx.get(), *sp_rhs.get());
            *sp_u_approx_tstart = csp_u_tstop_approx;
            if (fstop_ != nullptr) {
                auto fstop = *(SP_GridFunction*)fstop_->value;
                //ug::VecAdd(1.0, *sp_u_approx_tstart->get(), 1.0, *fstop.get());
                ug::VecAdd(*sp_u_approx_tstart->get(),*sp_u_approx_tstart->get(), *fstop.get());
            }
            return 0;
        };

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& status) override {
            this->Step(r_, u_, nullptr, status);
            this->Sum(1, u_, -1, r_);
            return 0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void print_settings() {}

        void setAdaptConv(bool conv) {}

        void setForceConv(bool force) {}

        void set_solver(SP_Solver solver) {
            this->linSolver = solver;
        }

        void set_approx_space(SP_ApproxSpace approx_space) {
            this->m_approx_space = approx_space;
        }

        //--------------------------------------------------------------------------------------------------------------

    };
}}*/
#endif