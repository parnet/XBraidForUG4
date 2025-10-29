#ifndef LIMEX_OBSERVER_HPP
#define LIMEX_OBSERVER_HPP

#include <vector>

#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"


namespace ug { namespace xbraid {
template <typename TDomain, typename TAlgebra>
class LimexObserver : public ITimeIntegratorObserver<TDomain,TAlgebra> {

    using T_GridFunction = GridFunction<TDomain, TAlgebra>;
    using SP_GridFunction = SmartPtr<T_GridFunction>;

    bool step_process(SP_GridFunction u, int step, number time, number dt){
        _step_u.push(u);
        _step_index.push_back(step);
        _step_time.push_back(time);
        _step_dt.push_back(dt);
        return true;
    };
    void clear();
    int size();

    std::vector<int> _step_index;
    std::vector<double> _step_time;
    std::vector<double> _step_dt;
    std::vector<SP_GridFunction> _step_u;
};


template<typename TDomain, typename TAlgebra>
int LimexObserver<TDomain,TAlgebra>::size(){
    return _step_index.size();
}
template<typename TDomain, typename TAlgebra>
void LimexObserver<TDomain,TAlgebra>::clear(){
    _step_index.clear();
    _step_time.clear();
    _step_dt.clear();
    _step_u.clear();
}


}}

#endif //LIMEX_OBSERVER_HPP
