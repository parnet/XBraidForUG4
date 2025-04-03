#ifndef UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_BDF_INTEGRATOR_H
#define UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_BDF_INTEGRATOR_H

#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BDF_Integrator : public ITimeIntegrator<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;
        using CP_GridFunction = ConstSmartPtr<T_GridFunction> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_VectorTimeSeries = VectorTimeSeries<typename TAlgebra::vector_type> ;
        using SP_VectorTimeSeries = SmartPtr<T_VectorTimeSeries> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_Operator = AssembledLinearOperator<TAlgebra> ;
        using SP_Operator = SmartPtr<T_Operator> ;

        using T_TimeDisc = BDF<TAlgebra> ;
        using SP_TimeDisc = SmartPtr<T_TimeDisc> ;

        //--------------------------------------------------------------------------------------------------------------

        BDF_Integrator() : ITimeIntegrator<TDomain, TAlgebra>() {}
        ~BDF_Integrator() override = default;

        //--------------------------------------------------------------------------------------------------------------


        void init(T_GridFunction const& u) override {};
        void prepare(T_GridFunction& u) override {};

        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) override {
            SP_GridFunction u0_nonconst = u0.cast_const();
            auto gridlevel = u0_nonconst->grid_level();
            if (!initialized_) {
                time_disc_ = make_sp(new T_TimeDisc(domain_disc_));
                operator_a_ = make_sp(new T_Operator(time_disc_, gridlevel));
                initialized_ = true;
            }
            bool success = false;

            auto solTimeSeries = make_sp(new T_VectorTimeSeries());
            solTimeSeries->push(u0_nonconst, t0);

            double current_dt = (t1 - t0) / this->order_;
            auto rhs = u0_nonconst->clone();
            double time = t0;

            time_disc_->set_order(1);
            for (int step = 0; step < this->order_; step++) {
                std::cout << "BDF Step  " << step + 1 << " / " << this->order_ << std::endl;
                time_disc_->prepare_step(solTimeSeries, current_dt);

                if (!assembled_ || fabs(assembled_dt_ - current_dt) > reassemble_threshold_) {
                    std::cout << "Assemble Operator for dt " << current_dt << " dt_prev=" << assembled_dt_
                        << std::endl;
                    time_disc_->assemble_jacobian(*operator_a_.get(), *u0_nonconst.get(), gridlevel);
                    assembled_ = true;
                    assembled_dt_ = current_dt;
                }
                time_disc_->assemble_rhs(*rhs.get(), gridlevel);
                linear_solver_->init(operator_a_, *u1.get());
                success = linear_solver_->apply(*u1.get(), *rhs.get());
                time += current_dt;
                if (step != this->order_ - 1) {
                    time_disc_->set_order(step + 1);
                    time = time_disc_->future_time();
                    solTimeSeries->push(u1->clone(), time);
                } else {
                    time = time_disc_->future_time();
                    solTimeSeries->push_discard_oldest(u1->clone(), time);
                }

            }
            return success;
        }

        //--------------------------------------------------------------------------------------------------------------

        void set_order(double order) {
            this->order_ = order;
        }

        void set_reassemble_threshold(double threshold) {
            this->reassemble_threshold_ = threshold;
        }
        void set_domain(SP_DomainDisc domain) {
            this->domain_disc_ = domain;
        }

        void set_solver(SP_Solver solver) {
            this->linear_solver_ = solver;
        }

        //--------------------------------------------------------------------------------------------------------------

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc domain_disc_;
        SP_Solver linear_solver_;
        SP_Operator operator_a_;
        SP_TimeDisc time_disc_;

        bool initialized_ = false;
        bool assembled_ = false;
        int order_ = 1;
        double reassemble_threshold_ = 1e-8;
        double assembled_dt_ = -1.0;

    };
}}

#endif