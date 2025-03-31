#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_VTK_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_VTK_OBSERVER_HPP


#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"
#include "lib_disc/io/vtkoutput.h"



namespace ug{ namespace xbraid {
    template<typename TDomain, typename TAlgebra>
    class VTK_Observer : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_VTKOutput = VTKOutput<TDomain::dim> ;
        using SP_VTKOutput = SmartPtr<T_VTKOutput> ;

        //--------------------------------------------------------------------------------------------------------------

        SP_VTKOutput m_out;
        const char *m_filename;

        //--------------------------------------------------------------------------------------------------------------

        VTK_Observer(SP_VTKOutput p_out, const char *filename) : ITimeIntegratorObserver<TDomain, TAlgebra>() {
            m_out = p_out;
            m_filename = filename;
        }

        ~VTK_Observer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            this->m_out->print(m_filename, *u, index, time);
            return true;
        };

        void write_time_pvd(SP_GridFunction u) {
            this->m_out->write_time_pvd(this->m_filename, *u);
        };

        void write_time_pvd_fn(const char *filename, SP_GridFunction u) {
            this->writeTimePVD(&u);
        };

        void set_filename(const char *filename) {
            this->m_filename = filename;
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif