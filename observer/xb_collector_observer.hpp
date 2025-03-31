#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_XB_COLLECTOR_PROCESS_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_XB_COLLECTOR_PROCESS_OBSERVER_HPP

#include "interface/observer_xbraid.hpp"


namespace ug{ namespace xbraid {
    template <typename TDomain, typename TAlgebra>
    class XBraidTimeIntegratorObserverCollector : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_TimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_TimeIntegratorObserver= SmartPtr<T_TimeIntegratorObserver> ;

        using T_XBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_XBraidTimeIntegratorObserver = SmartPtr<T_XBraidTimeIntegratorObserver> ;

        //--------------------------------------------------------------------------------------------------------------

        std::vector<SP_TimeIntegratorObserver> lst;
        std::vector<SP_XBraidTimeIntegratorObserver> xb_lst;

        //--------------------------------------------------------------------------------------------------------------

        XBraidTimeIntegratorObserverCollector() = default;

        ~XBraidTimeIntegratorObserverCollector() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void attach_observer(SP_XBraidTimeIntegratorObserver observer) {
            this->xb_lst.emplace_back(observer);
        }

        void attach_common_observer(SP_TimeIntegratorObserver observer) {
            this->lst.emplace_back(observer);
        }


        bool step_process(SP_GridFunction u, int index, number time, number dt) override {
            for (size_t k = 0; k < this->lst.size(); k++) {
                this->lst[k]->step_process(u, index, time,dt);
            }

            for (size_t k = 0; k < this->xb_lst.size(); k++) {
                this->xb_lst[k]->step_process(u, index, time,dt);
            }

            return true;
        };

        bool step_process(SP_GridFunction u, int index, number time, number dt , int iteration, int level) override {
            for (size_t k = 0; k < this->xb_lst.size(); k++) {
                this->xb_lst[k]->step_process(u, index, time, dt, iteration, level);
            }
            return true;
        };

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif