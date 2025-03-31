//todo check
#ifndef UGPLUGIN_XBRAIDFORUG4_FACTORY_SIMPLE_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_FACTORY_SIMPLE_INTEGRATOR_FACTORY_HPP

//2025-03 #include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"

#include "Limex/time_disc/linear_implicit_timestep.h"
#include "interface/integrator_factory.hpp"

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class SimpleIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_VectorType = typename TAlgebra::vector_type ;
        using T_ErrorEstimator = ISubDiagErrorEst<T_VectorType> ;
        using SP_ErrorEstimator = SmartPtr<T_ErrorEstimator> ;

        using T_IntegratorType = INonlinearTimeIntegrator<TDomain, TAlgebra> ;
        using T_SolverType = typename T_IntegratorType::solver_type;
        using SP_Solver = SmartPtr<T_SolverType> ;

        using T_ITimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_ITimeIntegrator = SmartPtr<T_ITimeIntegrator> ;


        using T_DomainDisc = IDomainDiscretization<TAlgebra> ;
        using SP_DomainDisc = SmartPtr<T_DomainDisc> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_DomainDisc m_domainDisc;
        SP_Solver m_solver;

        double m_dt_min = -1.0;
        double m_dt_max = -1.0;

        //--------------------------------------------------------------------------------------------------------------

        SimpleIntegratorFactory() : IntegratorFactory<TDomain, TAlgebra>() {        };

        ~SimpleIntegratorFactory() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_domain_disc(SP_DomainDisc dom) {
            this->m_domainDisc = dom;
        }

        void set_solver(SP_Solver solver) {
            this->m_solver = solver;
        }

        void set_dt_min(double dt_min) {
            this->m_dt_min = dt_min;
        }

        void set_dt_max(double dt_max) {
            this->m_dt_max = dt_max;
        }

        //--------------------------------------------------------------------------------------------------------------

        SP_ITimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) override {
            return create_time_integrator(current_dt, done);
        }

        SP_ITimeIntegrator create_time_integrator(double current_dt, bool done) override {
            auto time_stepper = make_sp(new LinearImplicitEuler<TAlgebra>(this->m_domainDisc));
            time_stepper->disable_matrix_cache();

            auto integrator = make_sp(new SimpleTimeIntegrator<TDomain, TAlgebra>(time_stepper));
            integrator->set_solver(this->m_solver);
            integrator->set_time_step(current_dt);
            integrator->set_dt_min(m_dt_min);
            integrator->set_dt_max(m_dt_max);
            return integrator;
        }
        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif