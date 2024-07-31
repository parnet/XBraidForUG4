#ifndef UG_PLUGIN_XBRAIDFORUG4_INTEGRATOR_FACTORY_H
#define UG_PLUGIN_XBRAIDFORUG4_INTEGRATOR_FACTORY_H

#include "../../../Limex/time_disc/time_integrator.hpp" // from Limex-Plugin


namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class IntegratorFactory {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef ITimeIntegrator<TDomain, TAlgebra> T_ITimeIntegrator;
        typedef SmartPtr<T_ITimeIntegrator> SP_ITimeIntegrator;

        //--------------------------------------------------------------------------------------------------------------

        IntegratorFactory() = default;

        virtual ~IntegratorFactory() = default;

        //--------------------------------------------------------------------------------------------------------------

        virtual SP_ITimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) = 0;

        virtual SP_ITimeIntegrator create_time_integrator(double current_dt, bool done) = 0;

        //--------------------------------------------------------------------------------------------------------------

    };

}}

#endif //UG_PLUGIN_XBRAIDFORUG4_INTEGRATOR_FACTORY_H
