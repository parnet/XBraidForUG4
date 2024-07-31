#ifndef UG_PLUGIN_XBRAIDFORUG4_VTK_PROCESS_OBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_VTK_PROCESS_OBSERVER_H

#include <map>
#include <tuple>

#include "lib_disc/io/vtkoutput.h"

#include "../interface/observer_xbraid.h"


namespace ug { namespace XBraidForUG4 {

    template<typename TDomain, typename TAlgebra>
    class VTK_ProcessObserver : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef VTKOutput<TDomain::dim> T_VTKOutput;
        typedef SmartPtr<T_VTKOutput> SP_VTKOutput;
        typedef std::tuple<int, int, int> TKey;

        //--------------------------------------------------------------------------------------------------------------

        SP_VTKOutput m_out;
        const char *m_filename;
        std::map<TKey, int> map;

        //--------------------------------------------------------------------------------------------------------------

        VTK_ProcessObserver(SP_VTKOutput p_out, const char *filename) : IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {
            m_out = p_out;
            m_filename = filename;
        }

        ~VTK_ProcessObserver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void write_time_pvd(SP_GridFunction u) { // todo neccessary? delte?
            this->m_out->write_time_pvd(this->m_filename, *u);
        };

        void write_time_pvd_fn(const char *filename, SP_GridFunction u) { // todo neccessary? delte?
            this->writeTimePVD(&u);
        };

        void set_filename(const char *filename) {
            this->m_filename = filename;
        }

        bool step_process(SP_GridFunction u, int index, number time, number dt) override {
            this->m_out->print(this->m_filename, *u, index, time);
            return true;
        };

        bool step_process(SP_GridFunction u, int index, double time, double dt, int iteration, int level) override {
            std::stringstream ss;
            int count = 0;
            auto tuple = std::make_tuple(index, iteration, level);
            auto it = this->map.find(tuple);
            if (it != this->map.end()) {
                count = it->second;
                count += 1;
                this->map[tuple] = count;
            } else {
                count = 0;
                this->map.emplace(tuple, 0);
            }
            ss << this->m_filename << "_k" << iteration << "_l" << level << "_c" << count;
            this->m_out->print(ss.str().c_str(), *u, index, time);
            return true;
        };



        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif //UG_PLUGIN_XBRAIDFORUG4_VTK_PROCESS_OBSERVER_H
