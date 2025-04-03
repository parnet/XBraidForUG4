#ifndef UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_THETA_INTEGRATOR_HPP
#define UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_THETA_INTEGRATOR_HPP

#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/solution_time_series.h"

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class ThetaIntegrator : public ITimeIntegrator<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;
        using CP_GridFunction = ConstSmartPtr<T_GridFunction> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SPDomainDisc = SmartPtr<T_DomainDisc> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_Operator = AssembledLinearOperator<TAlgebra> ;
        using SP_Operator = SmartPtr<T_Operator> ;

        using T_TimeDisc = ThetaTimeStep<TAlgebra> ;
        using SP_TimeDisc = SmartPtr<T_TimeDisc> ;

        using T_VectorTimeSeries = VectorTimeSeries<typename TAlgebra::vector_type> ;
        using SP_VectorTimeSeries = SmartPtr<T_VectorTimeSeries> ;



        //--------------------------------------------------------------------------------------------------------------

        ThetaIntegrator() : ITimeIntegrator<TDomain, TAlgebra>() {}
        ~ThetaIntegrator() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void init(T_GridFunction const& u) override {};

        void prepare(T_GridFunction& u) override {};


        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) override {

            auto solTimeSeries = make_sp(new T_VectorTimeSeries());
            SP_GridFunction u0_nonconst = u0.cast_const();

            auto gridlevel = u0_nonconst->grid_level();
            solTimeSeries->push(u0_nonconst, t0);
            double current_dt = t1 - t0;

            if (!initialized_) {
                time_disc_ = make_sp(new T_TimeDisc(domain_disc_));
                time_disc_->set_theta(theta);
                operator_a_ = make_sp(new T_Operator(time_disc_, gridlevel));
                initialized_ = true;
            }

            time_disc_->prepare_step(solTimeSeries, current_dt);
            // ø time_disc_->adjust_solution(*u0_nonconst.get(),gridlevel);

            auto rhs = u0_nonconst->clone();
            if (!assembled_ || fabs(assembled_dt_ - current_dt) > reassemble_threshold_) {
                std::cout << "Assemble Operator for dt " << current_dt << " dt_prev=" << assembled_dt_ << std::endl;
                time_disc_->assemble_jacobian(*operator_a_.get(), *u0_nonconst.get(), gridlevel);
                assembled_dt_ = current_dt;
                assembled_ = true;
            }
            time_disc_->assemble_rhs(*rhs.get(), gridlevel);
            // ø if (f != SPNULL) {
            // ø     (*rhs) += (*f);
            // ø  }

            linear_solver_->init(operator_a_, *u1.get());

            bool success = linear_solver_->apply(*u1.get(), *rhs.get());
            // ø solTimeSeries->push(u1,t1);
            // ø time_disc_->finish_step(solTimeSeries);
            return success;
        };

        //--------------------------------------------------------------------------------------------------------------

        void set_theta(double p_theta) {
            this->theta = p_theta;
        }
        void set_reassemble_threshold(double threshold) {
            this->reassemble_threshold = threshold;
        }
        void set_domain(SPDomainDisc domain) {
            this->m_domain_disc = domain;
        }

        void set_solver(SP_Solver solver) {
            this->m_linear_solver = solver;
        }
        //--------------------------------------------------------------------------------------------------------------

        SPDomainDisc domain_disc_;
        SP_Solver linear_solver_;

        SP_Operator operator_a_;
        SP_TimeDisc time_disc_;

        bool initialized_ = false;
        bool assembled_ = false;

        double theta = 1;

        double reassemble_threshold_ = 1e-8;
        double assembled_dt_ = -1.0;
        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif