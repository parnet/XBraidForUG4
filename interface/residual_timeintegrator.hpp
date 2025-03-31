#ifndef UGPLUGIN_XBRAIDFORUG4_INTERFACE_RESIDUAL_TIMEINTEGRATOR_HPP
#define UGPLUGIN_XBRAIDFORUG4_INTERFACE_RESIDUAL_TIMEINTEGRATOR_HPP


//2025-03 #include "../../../Limex/time_disc/interface/time_integrator.h"
#include "Limex/time_disc/time_integrator.hpp"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class IResidualTimeIntegrator : public ITimeIntegrator<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;
        using CP_GridFunction = ConstSmartPtr<T_GridFunction> ;

        //--------------------------------------------------------------------------------------------------------------

    protected:
        IResidualTimeIntegrator() : ITimeIntegrator<TDomain, TAlgebra>() {};
    public:
        ~IResidualTimeIntegrator() override = default;

        //--------------------------------------------------------------------------------------------------------------

        // virtual bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) = 0; derived from interface

        using ITimeIntegrator<TDomain,TAlgebra>::apply;

        virtual bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0, SP_GridFunction f) = 0;

        virtual SP_GridFunction defect(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) = 0;

        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif