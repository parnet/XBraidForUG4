#ifndef UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_THETA_CONSTSTEP_INTEGRATOR_NL_HPP
#define UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_THETA_CONSTSTEP_INTEGRATOR_NL_HPP


#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton.h"

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class ThetaConstStepIntegratorNL : public ITimeIntegrator<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------
        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;
        using CP_GridFunction = ConstSmartPtr<T_GridFunction> ;

        using T_DomainDisc =IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc =  SmartPtr<T_DomainDisc> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_NLSolver = NewtonSolver<TAlgebra> ;
        using SP_NLSolver = SmartPtr<T_NLSolver> ;

        using T_Operator = AssembledOperator<TAlgebra> ;
        using SP_Operator = SmartPtr<T_Operator> ;

        using T_TimeDisc = ThetaTimeStep<TAlgebra> ;
        using SP_TimeDisc = SmartPtr<T_TimeDisc> ;

        using T_VectorTimeSeries = VectorTimeSeries<typename TAlgebra::vector_type> ;
        using SP_VectorTimeSeries = SmartPtr<T_VectorTimeSeries> ;



        //--------------------------------------------------------------------------------------------------------------

        ThetaConstStepIntegratorNL() : ITimeIntegrator<TDomain, TAlgebra>() {}
        ~ThetaConstStepIntegratorNL() override = default;

        //--------------------------------------------------------------------------------------------------------------



        void init(T_GridFunction const& u) override {};

        void prepare(T_GridFunction& u) override {};

        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0,
                   number t0) override {
            SP_GridFunction u0_nonconst = u0.cast_const();
            auto gridlevel = u0_nonconst->grid_level();
            if (!initialized_) {
                time_disc_ = make_sp(new ThetaTimeStep<TAlgebra>(domain_disc_));
                time_disc_->set_theta(theta_);
                operator_a_ = make_sp(new AssembledOperator<TAlgebra>(time_disc_, gridlevel));
                initialized_ = true;
            }
            bool success = false;

            auto solTimeSeries = make_sp(new T_VectorTimeSeries());
            solTimeSeries->push(u0_nonconst, t0);
            non_linear_solver_->init(operator_a_);

            double current_dt = (t1 - t0) / this->num_steps_;
            auto ux = u0->clone();
            auto defaultLineSearch = non_linear_solver_->line_search();
            double time = t0;
            for (int step = 0; step < this->num_steps_; step++) {
                non_linear_solver_->set_line_search(defaultLineSearch);
                time_disc_->set_stage(1);
                time_disc_->prepare_step(solTimeSeries, current_dt);
                non_linear_solver_->prepare(*ux.get());
                success = non_linear_solver_->apply(*ux.get());
                if (step != this->num_steps_ - 1) {
                    time = time_disc_->future_time();
                    solTimeSeries->push_discard_oldest(ux->clone(), time);
                }
            }
            *u1 = *ux;
            return success;
        };

        //--------------------------------------------------------------------------------------------------------------

        void set_theta(double theta) {
            this->theta_ = theta;
        }

        void set_num_steps(double steps) {
            this->num_steps_ = steps;
        }

        void set_domain(SP_DomainDisc domain) {
            this->domain_disc_ = domain;
        }

        void set_solver(SP_NLSolver solver) {
            this->non_linear_solver_ = solver;
        }

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc domain_disc_;
        SP_NLSolver non_linear_solver_;
        SP_Operator operator_a_;
        SP_TimeDisc time_disc_;

        bool initialized_ = false;
        int num_steps_ = 1;
        double theta_ = 1;
        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif