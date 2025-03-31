//todo residual methods for step and defect could be implemented.

// check if the u_k are adequade for this method
// or timestep until the defect estimation
#ifndef UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_THETA_CONSTSTEP_INTEGRATOR_HPP
#define UGPLUGIN_XBRAIDFORUG4_INTEGRATOR_THETA_CONSTSTEP_INTEGRATOR_HPP

#include "lib_algebra/operator/linear_solver/linear_solver.h"

//2025-03 #include "lib_algebra/parallelization/parallel_vector_impl.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"

//2025-03 #include "math/vector.hpp"

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class ThetaConstStepIntegrator : public ITimeIntegrator<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;
        using CP_GridFunction = ConstSmartPtr<T_GridFunction> ;

        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        using T_Solver = LinearSolver<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_Solver> ;

        using T_Operator = AssembledLinearOperator<TAlgebra> ;
        using SP_Operator = SmartPtr<T_Operator> ;

        using T_TimeDisc = ThetaTimeStep<TAlgebra> ;
        using SP_TimeDisc = SmartPtr<T_TimeDisc> ;

        using T_VectorTimeSeries = VectorTimeSeries<typename TAlgebra::vector_type> ;
        using SP_VectorTimeSeries = SmartPtr<T_VectorTimeSeries> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc domain_disc_;
        SP_Solver linear_solver_;

        SP_Operator operator_a_;
        SP_TimeDisc time_disc_;

        bool initialized_ = false;
        bool assembled_ = false;

        int num_steps_ = 1;
        double theta_ = 1;

        double reassemble_threshold_ = 1e-8;
        double assembled_dt_ = -1.0;

        //--------------------------------------------------------------------------------------------------------------

        ThetaConstStepIntegrator() : ITimeIntegrator<TDomain, TAlgebra>() { }
        ~ThetaConstStepIntegrator() override = default;

        //--------------------------------------------------------------------------------------------------------------


        void init(T_GridFunction const& u) override {};

        void prepare(T_GridFunction& u) override {};

        bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) override {

            bool success = false;

            auto solTimeSeries = make_sp(new T_VectorTimeSeries());
            SP_GridFunction u0_nonconst = u0.cast_const();

            auto gridlevel = u0_nonconst->grid_level();
            solTimeSeries->push(u0_nonconst, t0);

            double current_dt = (t1 - t0) / this->num_steps;

            if (!initialized_) {
                time_disc_ = make_sp(new T_TimeDisc(domain_disc_));
                time_disc_->set_theta(theta_);
                operator_a_ = make_sp(new T_Operator(time_disc_, gridlevel));
                initialized_ = true;
            }

            auto rhs = u0_nonconst->clone();
            auto ux = u0->clone();

            double time = t0;
            for (int step = 0; step < this->num_steps; step++) {
                std::cout << "Theta step " << step + 1 << " / " << this->num_steps << std::endl;
                time_disc_->prepare_step(solTimeSeries, current_dt);
                if (!assembled_ || fabs(assembled_dt_ - current_dt) > reassemble_threshold_) {
                    std::cout << "Assemble Operator for dt " << current_dt << " dt_prev=" << assembled_dt_ << std::endl;
                    time_disc_->assemble_jacobian(*operator_a_, *u0_nonconst, gridlevel);
                    assembled_dt_ = current_dt;
                    assembled_ = true;
                }
                time_disc_->assemble_rhs(*rhs, gridlevel);

                //2025-03 double weight = 1.0 / (num_steps_ - step);
                //2025-03 VecScaleAdd(*ux, 1.0 - weight, *ux, weight, *u1);
                linear_solver_->init(operator_a_, *ux);
                success = linear_solver_->apply(*ux, *rhs);

                if (step != this->num_steps - 1) {
                    time += current_dt;
                    solTimeSeries->push(ux->clone(), time);
                }
            }
            *u1 = *ux;
            return success;
        };
        //--------------------------------------------------------------------------------------------------------------
        void set_theta(double p_theta) {
            this->theta = p_theta;
        }

        void set_num_steps(double steps) {
            this->num_steps = steps;
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
#endif