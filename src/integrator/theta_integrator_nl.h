#ifndef UG_PLUGIN_XBRAIDFORUG4_NL_THETA_INTEGRATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_NL_THETA_INTEGRATOR_H

#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"

namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class ThetaIntegratorNL : public ITimeIntegrator<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;
        typedef ConstSmartPtr<T_GridFunction> CP_GridFunction;

        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;

        typedef NewtonSolver<TAlgebra> T_NLSolver;
        typedef SmartPtr<T_NLSolver> SP_NLSolver;

        typedef AssembledOperator<TAlgebra> T_Operator;
        typedef SmartPtr<T_Operator> SP_Operator;

        typedef ThetaTimeStep<TAlgebra> T_TimeDisc;
        typedef SmartPtr<T_TimeDisc> SP_TimeDisc;

        typedef VectorTimeSeries<typename TAlgebra::vector_type> T_VectorTimeSeries;
        typedef SmartPtr<T_VectorTimeSeries> SP_VectorTimeSeries;

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc m_domain_disc;
        SP_NLSolver m_non_linear_solver;
        SP_Operator Operator_A;
        SP_TimeDisc time_disc;

        bool initialized = false;
        double theta = 1;

        //--------------------------------------------------------------------------------------------------------------

        ThetaIntegratorNL() : ITimeIntegrator<TDomain, TAlgebra>() {}
        ~ThetaIntegratorNL() override = default;

        //--------------------------------------------------------------------------------------------------------------


        void init(T_GridFunction const& u) override {};

        void prepare(T_GridFunction& u) override {};

        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) override {
            SP_GridFunction u0_nonconst = u0.cast_const()->clone();
            auto gridlevel = u0_nonconst->grid_level();
            if (!initialized) {
                time_disc = make_sp(new ThetaTimeStep<TAlgebra>(m_domain_disc));
                time_disc->set_theta(theta);
                Operator_A = make_sp(new AssembledOperator<TAlgebra>(time_disc, gridlevel));
                initialized = true;
                m_non_linear_solver->init(Operator_A);
                auto defaultLineSearch = m_non_linear_solver->line_search();
                m_non_linear_solver->set_line_search(defaultLineSearch);
            }
            bool success = false;
            auto solTimeSeries = make_sp(new T_VectorTimeSeries());
            solTimeSeries->push(u0_nonconst, t0);
            double current_dt = (t1 - t0);
            auto ux = u1->clone();
            time_disc->set_stage(1);
            time_disc->prepare_step(solTimeSeries, current_dt);
            auto result = m_non_linear_solver->prepare(*ux.get());
            std::cout << "prepare " << result << std::endl << std::flush;
            success = m_non_linear_solver->apply(*ux.get());
            *u1 = *ux;
            return success;
        };

        //--------------------------------------------------------------------------------------------------------------

        void set_theta(double p_theta) {
            this->theta = p_theta;
        }

        void set_domain(SP_DomainDisc domain) {
            this->m_domain_disc = domain;
        }

        void set_solver(SP_NLSolver solver) {
            this->m_non_linear_solver = solver;
        }
        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif //UG_PLUGIN_XBRAIDFORUG4_NL_THETA_INTEGRATOR_H
