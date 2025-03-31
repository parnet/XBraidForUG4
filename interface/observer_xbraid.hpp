#ifndef UGPLUGIN_XBRAIDFORUG4_INTERFACE_XBRAID_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_INTERFACE_XBRAID_OBSERVER_HPP

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"

namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class IXBraidTimeIntegratorObserver  : public ITimeIntegratorObserver<TDomain, TAlgebra>
                     /*public ITimeIntegratorStageObserver_start<TDomain, TAlgebra>,
                     public ITimeIntegratorStageObserver_finalize<TDomain, TAlgebra>,
                     public ITimeIntegratorStageObserver_end<TDomain, TAlgebra>*/
    {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        //--------------------------------------------------------------------------------------------------------------
    protected:
        IXBraidTimeIntegratorObserver() = default;
    public:
        ~IXBraidTimeIntegratorObserver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * defines a callback which is called in specified states by the method where the observer is attached
         * see also ITimeIntegratorObserver documentation
         * @param u the gridfunction solution from the method
         * @param index the index e.g. the number of timesteps that have already been processed
         * @param time the timepoint at which this solution is evaluated
         * @param dt  the time stepsize which is used
         * @return
         */
        bool step_process(SP_GridFunction u, int index, number time, number dt) override = 0;

        /**
         * defines a callback which is called in specified states by the method where the observer is attached
         * in addition to the normal observer method this method can be called by additional states and provides
         * information about the level and the current iteration of the method
         * @param u the gridfunction solution from the method
         * @param index the index e.g. the number of timesteps that have already been processed
         * @param time the timepoint at which this solution is evaluated
         * @param dt the time stepsize which is used
         * @param iteration the current iteration of the method
         * @param level the grid level
         * @return
         */
        virtual bool step_process(SP_GridFunction u, int index, number time, number dt, int iteration, int level) = 0;

        /**
         * [[Depcrated]] A helper method that is registered to be used in Lua to access the other
         * @param u the gridfunction solution from the method
         * @param index the index e.g. the number of timesteps that have already been processed
         * @param time the timepoint at which this solution is evaluated
         * @param dt the time stepsize which is used
         * @return
         */
        bool write(SP_GridFunction u, int index, number time, number dt) {
            return this->step_process(u,index,time,dt);
        }


        //--------------------------------------------------------------------------------------------------------------
    };

 }}

#endif