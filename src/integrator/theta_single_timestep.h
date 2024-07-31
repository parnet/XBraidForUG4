#ifndef UG_PLUGIN_XBRAIDFORUG4_THETA_STEPPER_H
#define UG_PLUGIN_XBRAIDFORUG4_THETA_STEPPER_H

#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "../interface/residual_timeintegrator.h"



namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class ThetaSingleTimeStep : public IResidualTimeIntegrator<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;
        typedef ConstSmartPtr<T_GridFunction> CP_GridFunction;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef LinearSolver<typename TAlgebra::vector_type> T_Solver;
        typedef SmartPtr<T_Solver> SP_Solver;

        typedef AssembledLinearOperator<TAlgebra> T_Operator;
        typedef SmartPtr<T_Operator> SP_Operator;

        typedef ThetaTimeStep<TAlgebra> T_TimeDisc;
        typedef SmartPtr<T_TimeDisc> SP_TimeDisc;

        typedef VectorTimeSeries<typename TAlgebra::vector_type> T_VectorTimeSeries;
        typedef SmartPtr<T_VectorTimeSeries> SP_VectorTimeSeries;

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc m_domain_disc;
        SP_Solver m_linear_solver;

        SP_Operator Operator_A;
        SP_TimeDisc time_disc;

        bool initialized = false;
        bool assembled = false;
        double theta = 1;
        double reassemble_threshold = 1e-8;
        double assembled_dt = -1.0;

        //--------------------------------------------------------------------------------------------------------------

        ThetaSingleTimeStep() : IResidualTimeIntegrator<TDomain, TAlgebra>() {}
        ~ThetaSingleTimeStep() override = default;

        //--------------------------------------------------------------------------------------------------------------


        void init(T_GridFunction const& u) override {};

        void prepare(T_GridFunction& u) override {};

        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) override {
            return this->apply(u1, t1, u0, t0, SPNULL);
        };

        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0, SP_GridFunction f) override {
            bool success = false;
            auto solTimeSeries = make_sp(new T_VectorTimeSeries());
            SmartPtr<T_GridFunction> u0_nonconst = u0.cast_const();
            auto gridlevel = u0_nonconst->grid_level();
            solTimeSeries->push(u0_nonconst, t0);
            double current_dt = t1 - t0;
            if (!initialized) {
                time_disc = make_sp(new T_TimeDisc(m_domain_disc));
                time_disc->set_theta(theta);
                Operator_A = make_sp(new T_Operator(time_disc, gridlevel));
                initialized = true;
            } // todo abstract and combine with defect method?
            time_disc->prepare_step(solTimeSeries, current_dt);
            auto rhs = u0_nonconst->clone();
            if (!assembled || fabs(assembled_dt - current_dt) > reassemble_threshold) {
                std::cout << "Assemble Operator for dt " << current_dt << " dt_prev=" << assembled_dt << std::endl;
                time_disc->assemble_jacobian(*Operator_A.get(), *u0_nonconst.get(), gridlevel);
                assembled_dt = current_dt;
                assembled = true;
            }
            time_disc->assemble_rhs(*rhs.get(), gridlevel);
            if (f != SPNULL) {
                (*rhs) += (*f); // todo einzige zeile die unterschiedlich ist, abstract !
            }
            m_linear_solver->init(Operator_A, *u1.get());
            success = m_linear_solver->apply(*u1.get(), *rhs.get());
            return success;
        };

        SP_GridFunction defect(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) override {
            bool success = false;
            auto solTimeSeries = make_sp(new ug::VectorTimeSeries<typename TAlgebra::vector_type>());
            SP_GridFunction u0_nonconst = u0.cast_const();
            auto gridlevel = u0_nonconst->grid_level();
            solTimeSeries->push(u0_nonconst, t0);
            double current_dt = t1 - t0;
            if (!initialized) {
                time_disc = make_sp(new ug::ThetaTimeStep<TAlgebra>(m_domain_disc));
                time_disc->set_theta(theta);
                Operator_A = make_sp(new ug::AssembledLinearOperator<TAlgebra>(time_disc, gridlevel));
                initialized = true;
            }
            time_disc->prepare_step(solTimeSeries, current_dt);
            auto rhs = u0_nonconst->clone();
            if (!assembled || fabs(assembled_dt - current_dt) > reassemble_threshold) {
                std::cout << "Assemble Operator for dt " << current_dt << " dt_prev=" << assembled_dt << std::endl;
                time_disc->assemble_jacobian(*Operator_A.get(), *u0_nonconst.get(), gridlevel);
                assembled_dt = current_dt;
                assembled = true;
            }
            time_disc->assemble_rhs(*rhs.get(), gridlevel);
            this->m_linear_solver->init(Operator_A, *u0);
            Operator_A->apply_sub(
                *rhs.get(), // f co domain function [in / out]
                *u1.get() // u domain function [in]
            ); // calculates r = r - A * u
            return rhs;
        };

        //--------------------------------------------------------------------------------------------------------------

        void set_theta(double p_theta) {
            this->theta = p_theta;
        }

        void set_reassemble_threshold(double threshold) {
            this->reassemble_threshold = threshold;
        }

        void set_domain(SP_DomainDisc domain) {
            this->m_domain_disc = domain;
        }

        void set_solver(SP_Solver solver) {
            this->m_linear_solver = solver;
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif //UG_PLUGIN_XBRAIDFORUG4_THETA_STEPPER_H
