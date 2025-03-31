#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_COLLECTOR_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_COLLECTOR_OBSERVER_HPP

#include <vector>

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class TimeIntegratorObserverCollector : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;


        using T_ITimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_ITimeIntegratorObserver = SmartPtr<T_ITimeIntegratorObserver> ;

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_ITimeIntegratorObserver> lst;

        //--------------------------------------------------------------------------------------------------------------

        TimeIntegratorObserverCollector() = default;
        ~TimeIntegratorObserverCollector() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void attach_observer(SP_ITimeIntegratorObserver observer) {
            this->lst.emplace_back(observer);
        }


        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            bool overall_result = true;
            for (size_t k = 0; k < this->lst.size(); k++) {
                bool this_result = this->lst[k]->step_process(u, index, time, dt);
                overall_result = overall_result && this_result;
            }
            return overall_result;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif