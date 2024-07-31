#ifndef UG_PLUGIN_XBRAIDFORUG4_MOCK_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_MOCK_OBSERVER_H


#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"


namespace ug { namespace XBraidForUG4 {

    template<typename TDomain, typename TAlgebra>
    class NoObserver : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        NoObserver() : ITimeIntegratorObserver<TDomain, TAlgebra>() {}

        ~NoObserver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            return true;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_MOCK_OBSERVER_H
