#ifndef UGPLUGIN_XBRAIDFORUG4_FACTORY_CONST_STEP_LINEAR_TIME_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_FACTORY_CONST_STEP_LINEAR_TIME_INTEGRATOR_FACTORY_HPP

//2025-03 #include "../../Limex/time_disc/const_step_linear_time_integrator.h"
#include "interface/integrator_factory.hpp"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class ConstStepLinearTimeIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------
        using T_SolverType = ILinearOperatorInverse<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_SolverType> ;

        using T_TimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_ConstStepLinearTimeIntegrator = ConstStepLinearTimeIntegrator<TDomain, TAlgebra> ;
        using SP_ConstStepLinearTimeIntegrator = SmartPtr<T_ConstStepLinearTimeIntegrator> ;

        using T_TimeDisc = ITimeDiscretization<TAlgebra> ;
        using SP_TimeDisc = SmartPtr<T_TimeDisc> ;



        //--------------------------------------------------------------------------------------------------------------

        ConstStepLinearTimeIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {};
        ~ConstStepLinearTimeIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_num_steps(int steps) {
            this->num_steps_ = steps;
        }

        void set_time_disc(SP_TimeDisc time_disc) {
            this->time_disc_ = time_disc;
        }

        void set_solver(SP_Solver solver) {
            this->solver_ = solver;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            return create_time_integrator(current_dt, done);
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            auto integrator = make_sp( new T_ConstStepLinearTimeIntegrator(time_disc_, solver_));
            integrator->set_num_steps(this->num_steps_);
            return integrator;
        }
        //--------------------------------------------------------------------------------------------------------------

        SP_Solver solver_;
        SP_TimeDisc time_disc_;
        int num_steps_ = 1;

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif