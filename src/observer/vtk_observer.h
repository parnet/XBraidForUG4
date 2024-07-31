#ifndef UG_PLUGIN_XBRAIDFORUG4_VTK_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_VTK_OBSERVER_H


#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"
#include "lib_disc/io/vtkoutput.h"



namespace ug { namespace XBraidForUG4 {
    template<typename TDomain, typename TAlgebra>
    class VTK_Observer : public ITimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef VTKOutput<TDomain::dim> T_VTKOutput;
        typedef SmartPtr<T_VTKOutput> SP_VTKOutput;

        //--------------------------------------------------------------------------------------------------------------

        SP_VTKOutput m_out;
        const char *m_filename;

        //--------------------------------------------------------------------------------------------------------------

        VTK_Observer(SP_VTKOutput p_out, const char *filename) : ITimeIntegratorObserver<TDomain, TAlgebra>() {
            m_out = p_out;
            m_filename = filename;
        }

        ~VTK_Observer() override {};

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
#endif //UG_PLUGIN_XBRAIDFORUG4_VTK_OBSERVER_H
