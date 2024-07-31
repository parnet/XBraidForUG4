//todo check
#ifndef UG_PLUGIN_XBRAIDFORUG4_SIMPLE_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_SIMPLE_INTEGRATOR_FACTORY_H

#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"

#include "../../../Limex/time_disc/linear_implicit_timestep.h"
#include "../../../Limex/time_disc/time_integrator.hpp"
#include "../../../Limex/time_disc/time_extrapolation.h"

#include "../interface/integrator_factory.h"


namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class SimpleIntegratorFactory : public IntegratorFactory<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef typename TAlgebra::vector_type T_VectorType;
        typedef ISubDiagErrorEst<T_VectorType> T_ErrorEstimator;
        typedef SmartPtr<T_ErrorEstimator> SP_ErrorEstimator;

        typedef INonlinearTimeIntegrator<TDomain, TAlgebra> T_IntegratorType;
        typedef typename T_IntegratorType::solver_type T_SolverType;
        typedef SmartPtr<T_SolverType> SP_Solver;

        typedef ITimeIntegrator<TDomain, TAlgebra> T_ITimeIntegrator;
        typedef SmartPtr<T_ITimeIntegrator> SP_ITimeIntegrator;


        typedef IDomainDiscretization<TAlgebra> T_DomainDisc;
        typedef SmartPtr<T_DomainDisc> SP_DomainDisc;




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
#endif //UG_PLUGIN_XBRAIDFORUG4_SIMPLE_INTEGRATOR_FACTORY_H
