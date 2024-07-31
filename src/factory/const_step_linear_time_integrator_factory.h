//todo check
#ifndef UG_PLUGIN_XBRAIDINTEGRATOR_CONST_STEP_LINEAR_TIME_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDINTEGRATOR_CONST_STEP_LINEAR_TIME_INTEGRATOR_FACTORY_H

#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"

#include "../../../Limex/time_disc/time_integrator.hpp"
#include "../../../Limex/time_disc/time_extrapolation.h"

#include "../interface/integrator_factory.h"



namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class ConstStepLinearTimeIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------
        typedef ILinearOperatorInverse<typename TAlgebra::vector_type> T_SolverType;
        typedef SmartPtr<T_SolverType> SP_Solver;

        typedef ITimeIntegrator<TDomain, TAlgebra> T_TimeIntegrator;
        typedef SmartPtr<T_TimeIntegrator> SP_TimeIntegrator;

        typedef ConstStepLinearTimeIntegrator<TDomain, TAlgebra> T_ConstStepLinearTimeIntegrator;
        typedef SmartPtr<T_ConstStepLinearTimeIntegrator> SP_ConstStepLinearTimeIntegrator;

        typedef ITimeDiscretization<TAlgebra> T_TimeDisc;
        typedef SmartPtr<T_TimeDisc> SP_TimeDisc;

        //--------------------------------------------------------------------------------------------------------------

        SP_Solver m_solver;
        SP_TimeDisc m_time_disc;
        int m_num_steps = 1;

        //--------------------------------------------------------------------------------------------------------------

        ConstStepLinearTimeIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {};
        ~ConstStepLinearTimeIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_num_steps(int steps) {
            this->m_num_steps = steps;
        }

        void set_time_disc(SP_TimeDisc time_disc) {
            this->m_time_disc = time_disc;
        }

        void set_solver(SP_Solver solver) {
            this->m_solver = solver;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            return create_time_integrator(current_dt, done);
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            auto integrator = make_sp( new T_ConstStepLinearTimeIntegrator(m_time_disc, m_solver));
            integrator->set_num_steps(this->m_num_steps);
            return integrator;
        }
        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDINTEGRATOR_CONST_STEP_LINEAR_TIME_INTEGRATOR_FACTORY_H
