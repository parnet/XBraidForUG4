#ifndef UG_PLUGIN_XBRAIDFORUG4_TIME_STEPPER_H
#define UG_PLUGIN_XBRAIDFORUG4_TIME_STEPPER_H

#include "../../../Limex/time_disc/time_integrator.hpp" // from Limex-Plugin

namespace ug { namespace XBraidForUG4 {
    template <typename TDomain, typename TAlgebra>
    class IResidualTimeIntegrator : public ITimeIntegrator<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;
        typedef ConstSmartPtr<T_GridFunction> CP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        IResidualTimeIntegrator() : ITimeIntegrator<TDomain, TAlgebra>() {};
        ~IResidualTimeIntegrator() override = default;

        //--------------------------------------------------------------------------------------------------------------

        virtual bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) = 0;

        virtual bool apply(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0, SP_GridFunction f) = 0;

        virtual SP_GridFunction defect(SP_GridFunction u1, number t1, CP_GridFunction u0, number t0) = 0;

        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif //UG_PLUGIN_XBRAIDFORUG4_TIME_STEPPER_H
