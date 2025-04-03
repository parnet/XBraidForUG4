#ifndef UGPLUGIN_XBRAIDFORUG4_OBSERVER_VTK_PROCESS_OBSERVER_HPP
#define UGPLUGIN_XBRAIDFORUG4_OBSERVER_VTK_PROCESS_OBSERVER_HPP

#include <map>
#include <tuple>

#include "lib_disc/io/vtkoutput.h"

#include "interface/observer_xbraid.hpp"


namespace ug{ namespace xbraid {

    template<typename TDomain, typename TAlgebra>
    class VTK_ProcessObserver : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_VTKOutput = VTKOutput<TDomain::dim> ;
        using SP_VTKOutput = SmartPtr<T_VTKOutput> ;
        using TKey = std::tuple<int, int, int> ;


        //--------------------------------------------------------------------------------------------------------------

        VTK_ProcessObserver(SP_VTKOutput p_out, const char *filename) : IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {
            out_ = p_out;
            filename_ = filename;
        }

        ~VTK_ProcessObserver() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void write_time_pvd(SP_GridFunction u) {
            this->m_out->write_time_pvd(this->m_filename, *u);
        };

        void write_time_pvd_fn(const char *filename, SP_GridFunction u) {
            this->writeTimePVD(&u);
        };

        void set_filename(const char *filename) {
            this->filename_ = filename;
        }

        bool step_process(SP_GridFunction u, int index, number time, number dt) override {
            this->out_->print(this->filename_, *u, index, time);
            return true;
        };

        bool step_process(SP_GridFunction u, int index, double time, double dt, int iteration, int level) override {
            std::stringstream ss;
            int count = 0;
            auto tuple = std::make_tuple(index, iteration, level);
            auto it = this->map_.find(tuple);
            if (it != this->map_.end()) {
                count = it->second;
                count += 1;
                this->map_[tuple] = count;
            } else {
                count = 0;
                this->map_.emplace(tuple, 0);
            }
            ss << this->filename_ << "_k" << iteration << "_l" << level << "_c" << count;
            this->out_->print(ss.str().c_str(), *u, index, time);
            return true;
        };


        //--------------------------------------------------------------------------------------------------------------

        SP_VTKOutput out_;
        const char *filename_;
        std::map<TKey, int> map_;

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif