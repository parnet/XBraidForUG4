#ifndef UG_PLUGIN_XBRAIDFORUG4_BDF_INTEGRATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_BDF_INTEGRATOR_H

#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"



namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BDF_Integrator : public ITimeIntegrator<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;
        typedef ConstSmartPtr<T_GridFunction> CP_GridFunction;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef VectorTimeSeries<typename TAlgebra::vector_type> T_VectorTimeSeries;
        typedef SmartPtr<T_VectorTimeSeries> SP_VectorTimeSeries;

        typedef LinearSolver<typename TAlgebra::vector_type> T_Solver;
        typedef SmartPtr<T_Solver> SP_Solver;

        typedef AssembledLinearOperator<TAlgebra> T_Operator;
        typedef SmartPtr<T_Operator> SP_Operator;

        typedef BDF<TAlgebra> T_TimeDisc;
        typedef SmartPtr<T_TimeDisc> SP_TimeDisc;

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc m_domain_disc;
        SP_Solver m_linear_solver;
        SP_Operator Operator_A;
        SP_TimeDisc time_disc;

        bool initialized = false;
        bool assembled = false;
        int order = 1;
        double reassemble_threshold = 1e-8;
        double assembled_dt = -1.0;

        //--------------------------------------------------------------------------------------------------------------

        BDF_Integrator() : ITimeIntegrator<TDomain, TAlgebra>() {}
        ~BDF_Integrator() override = default;

        //--------------------------------------------------------------------------------------------------------------


        void init(T_GridFunction const& u) override {};
        void prepare(T_GridFunction& u) override {};

        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) override {
            SP_GridFunction u0_nonconst = u0.cast_const();
            auto gridlevel = u0_nonconst->grid_level();
            if (!initialized) {
                time_disc = make_sp(new T_TimeDisc(m_domain_disc));
                Operator_A = make_sp(new T_Operator(time_disc, gridlevel));
                initialized = true;
            }
            bool success = false;
            auto solTimeSeries = make_sp(new T_VectorTimeSeries());
            solTimeSeries->push(u0_nonconst, t0);
            double current_dt = (t1 - t0) / this->order;
            auto rhs = u0_nonconst->clone();
            double time = t0;
            time_disc->set_order(1);
            for (int step = 0; step < this->order; step++) {
                std::cout << "BDF Step  " << step + 1 << " / " << this->order << std::endl;
                time_disc->prepare_step(solTimeSeries, current_dt);
                if (!assembled || fabs(assembled_dt - current_dt) > reassemble_threshold) {
                    std::cout << "Assemble Operator for dt " << current_dt << " dt_prev=" << assembled_dt
                        << std::endl;
                    time_disc->assemble_jacobian(*Operator_A.get(), *u0_nonconst.get(), gridlevel);
                    assembled = true;
                    assembled_dt = current_dt;
                }
                time_disc->assemble_rhs(*rhs.get(), gridlevel);
                m_linear_solver->init(Operator_A, *u1.get());
                success = m_linear_solver->apply(*u1.get(), *rhs.get());
                time += current_dt;
                time_disc->set_order(step + 1);
                solTimeSeries->push(u1->clone(), time);
            }
            return success;
        }

        //--------------------------------------------------------------------------------------------------------------

        void set_order(double p_order) {
            this->order = p_order;
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

#endif //UG_PLUGIN_XBRAIDFORUG4_BDF_INTEGRATOR_H
