//todo check || composite conv check => remove
#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAID_RESIDUAL_STEPPER_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAID_RESIDUAL_STEPPER_H

#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/operator/composite_conv_check.h"
#include "gridfunction_base.h"
#include "../../../Limex/time_disc/linear_implicit_timestep.h"
#include "../../../Limex/time_disc/time_extrapolation.h"
#include "../core/space_time_communicator.h"
#include "../util/parallel_logger.h"


namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BraidResidualStepper : public BraidGridFunctionBase<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef typename TAlgebra::vector_type T_Vector;

        typedef ISubDiagErrorEst<T_Vector> T_ErrorEstimatorType;
        typedef SmartPtr<T_ErrorEstimatorType> SP_ErrorEstimatorType;

        typedef AitkenNevilleTimex<T_Vector> T_TimexType;
        typedef SmartPtr<T_TimexType> SP_TimexType;

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef ThetaTimeStep<TAlgebra> T_TimeStep;
        typedef SmartPtr<T_TimeStep> SP_TimeStep;

        typedef INonlinearTimeIntegrator<TDomain, TAlgebra> T_INonlinearTimeIntegrator;
        typedef typename T_INonlinearTimeIntegrator::solver_type T_SolverType;
        typedef SmartPtr<T_SolverType> SP_SolverType;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef SimpleTimeIntegrator<TDomain, TAlgebra> T_TimeIntegratorType;
        typedef SmartPtr<T_TimeIntegratorType> SP_TimeIntegratorType;

        typedef LinearImplicitEuler<TAlgebra> T_TimeStepType;
        typedef SmartPtr<T_TimeStepType> SP_TimeStepType;

        typedef StdConvCheck<typename TAlgebra::vector_type> T_Conv;
        typedef SmartPtr<T_Conv> SP_Conv;

        typedef SmartPtr<SpaceTimeCommunicator> SP_SpaceTimeCommunicator;

        typedef SmartPtr<ParallelLogger> SP_Paralog;

        typedef LinearSolver<typename TAlgebra::vector_type> T_Solver;
        typedef SmartPtr<T_Solver> SP_Solver;

        typedef ApproximationSpace<TDomain> T_ApproxSpace;
        typedef SmartPtr<T_ApproxSpace> SP_ApproxSpace;

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


        int Step(braid_Vector u_, braid_Vector ustop_, braid_Vector fstop_, BraidStepStatus& pstatus) override {
            //print_status(this->m_log->o,pstatus);
            int level;
            int idone;
            int tindex;
            int iteration;
            double t_start;
            double t_stop;
            pstatus.GetLevel(&level);
            pstatus.GetTstartTstop(&t_start, &t_stop);
            pstatus.GetDone(&idone);
            pstatus.GetTIndex(&tindex);
            pstatus.GetIter(&iteration);
            double current_dt = t_stop - t_start;

            auto* sp_u_approx_tstart = (SP_GridFunction*)u_->value;
            auto* sp_u_tstop_approx = (SP_GridFunction*)ustop_->value;
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
                ::VecAdd(1.0, *sp_u_approx_tstart->get(), 1.0, *fstop.get());
            }
            return 0;
        };

        int Residual(braid_Vector u_, braid_Vector r_, BraidStepStatus& pstatus) override {
            this->Step(r_, u_, nullptr, pstatus);
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
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_BRAID_RESIDUAL_STEPPER_H
