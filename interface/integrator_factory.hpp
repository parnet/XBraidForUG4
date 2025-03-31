#ifndef UGPLUGIN_XBRAIDFORUG4_INTERFACE_INTEGRATOR_FACTORY_HPP
#define UGPLUGIN_XBRAIDFORUG4_INTERFACE_INTEGRATOR_FACTORY_HPP

#include "Limex/time_disc/time_integrator.hpp" // from Limex-Plugin





namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class IntegratorFactory {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_ITimeIntegrator = ITimeIntegrator<TDomain, TAlgebra> ;
        using SP_ITimeIntegrator = SmartPtr<T_ITimeIntegrator> ;

        //--------------------------------------------------------------------------------------------------------------
    protected:
        IntegratorFactory() = default;
    public:
        virtual ~IntegratorFactory() = default;

        //--------------------------------------------------------------------------------------------------------------

        virtual SP_ITimeIntegrator create_level_time_integrator(double current_dt, bool done, int level) = 0;

        virtual SP_ITimeIntegrator create_time_integrator(double current_dt, bool done) = 0;

        //--------------------------------------------------------------------------------------------------------------

    };

}}

#endif