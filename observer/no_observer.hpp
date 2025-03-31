#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_MOCK_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_MOCK_OBSERVER_HPP


#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"


/* 2025-03-25

namespace ug{ namespace xbraid {

    template<typename TDomain, typename TAlgebra>
    class NoObserver : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

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
*/
#endif