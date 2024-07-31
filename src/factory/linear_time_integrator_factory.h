#ifndef UG_PLUGIN_XBRAIDINTEGRATOR_LINEAR_TIME_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDINTEGRATOR_LINEAR_TIME_INTEGRATOR_FACTORY_H


#include "../../../Limex/time_disc/time_integrator.hpp"
#include "../../../Limex/time_disc/time_extrapolation.h"

#include "../interface/integrator_factory.h"

namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class LinearTimeIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef ILinearOperatorInverse<typename TAlgebra::vector_type> T_SolverType;
        typedef SmartPtr<T_SolverType> SP_Solver;

        typedef ITimeIntegrator<TDomain, TAlgebra> T_TimeIntegrator;
        typedef SmartPtr<T_TimeIntegrator> SP_TimeIntegrator;

        typedef ITimeDiscretization<TAlgebra> T_TimeDisc;
        typedef SmartPtr<T_TimeDisc> SP_TimeDisc;


        //--------------------------------------------------------------------------------------------------------------


        SP_Solver m_solver;
        SP_TimeDisc m_time_disc;

        //--------------------------------------------------------------------------------------------------------------

        LinearTimeIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {};
        ~LinearTimeIntegratorFactory() = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_time_disc(SP_TimeDisc tdisc) {
            this->m_time_disc = tdisc;
        }

        void set_solver(SP_Solver solver) {
            this->m_solver = solver;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            return create_time_integrator(current_dt, done);
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            auto integrator = make_sp(
                new LinearTimeIntegrator<TDomain, TAlgebra>(this->m_time_disc, this->m_solver));
            integrator->set_time_step(current_dt);

            return integrator;
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif //UG_PLUGIN_XBRAIDINTEGRATOR_LINEAR_TIME_INTEGRATOR_FACTORY_H
