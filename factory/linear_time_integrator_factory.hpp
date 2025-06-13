#ifndef UGPLUGIN_XBRAIDFORUG4_FACTORY_LINEAR_TIME_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_FACTORY_LINEAR_TIME_INTEGRATOR_FACTORY_HPP


//#include "../../../Limex/time_disc/linear_integrator/linear_time_integrator.h"

#include "interface/integrator_factory.hpp"

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class LinearTimeIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_SolverType = ILinearOperatorInverse<typename TAlgebra::vector_type> ;
        using SP_Solver = SmartPtr<T_SolverType> ;

        using T_TimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_TimeIntegrator = SmartPtr<T_TimeIntegrator> ;

        using T_TimeDisc = ITimeDiscretization<TAlgebra> ;
        using SP_TimeDisc = SmartPtr<T_TimeDisc> ;




        //--------------------------------------------------------------------------------------------------------------

        LinearTimeIntegratorFactory() = default;
        ~LinearTimeIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_time_disc(SP_TimeDisc tdisc) {
            this->time_disc_ = tdisc;
        }

        void set_solver(SP_Solver solver) {
            this->solver_ = solver;
        }

        SP_TimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            return create_time_integrator(current_dt, done);
        }

        SP_TimeIntegrator create_time_integrator(double current_dt, bool done) override {
            auto integrator = make_sp(
                new LinearTimeIntegrator<TDomain, TAlgebra>(this->time_disc_, this->solver_));
            integrator->set_time_step(current_dt);

            return integrator;
        }

        //--------------------------------------------------------------------------------------------------------------

        SP_Solver solver_;
        SP_TimeDisc time_disc_;
        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif