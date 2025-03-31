#ifndef UGPLUGIN_XBRAIDFORUG4_BRAID_HEAT_CHECK_HPP
#define UGPLUGIN_XBRAIDFORUG4_BRAID_HEAT_CHECK_HPP



#include "interface/observer_xbraid.hpp"
#include "observer/io_process_observer.hpp"
#include "observer/vtk_process_observer.hpp"
#include "util/parallel_logger.hpp"
//2025-03  #include "util/parallel_io_gridfunction.h"
//2025-03  #include "braid_biot_estimator.h"
#include "biot_error_data.hpp"
//2025-03 #include <stdint.h>
#include <lib_disc/function_spaces/interpolate.h>

namespace ug { namespace xbraid {namespace poro {

    template<typename TDomain, typename TAlgebra>
    class BraidHeatCheck : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_VTKOutput = VTKOutput<TDomain::dim> ;
        using SP_VTKOutput = SmartPtr<T_VTKOutput> ;

        using T_VTK_ProcessObserver = VTK_ProcessObserver<TDomain, TAlgebra> ;
        using SP_VTK_ProcessObserver = SmartPtr<T_VTK_ProcessObserver> ;

        using T_IO_ProcessObserver = IO_ProcessObserver<TDomain, TAlgebra> ;
        using SP_IO_ProcessOobserver = SmartPtr<T_IO_ProcessObserver> ;

        using T_ParallelLogger = ParallelLogger ;
        using SP_ParallelLogger = SmartPtr<T_ParallelLogger> ;
        using T_Key= std::tuple<int, int, int> ;

        //--------------------------------------------------------------------------------------------------------------

        BiotErrorData<TDomain, TAlgebra> err_u;
        BiotErrorData<TDomain, TAlgebra> err_sol;
        BiotErrorData<TDomain, TAlgebra> err_udiffsol;

        SP_VTK_ProcessObserver m_out_solution;
        SP_VTK_ProcessObserver m_out_diff;

        SP_IO_ProcessOobserver m_ioout_solution;
        SP_IO_ProcessOobserver m_ioout_diff;

        SP_ParallelLogger m_log;

        int num_ref = 3;
        int max_index = 512;
        int max_index_precomputed = 512;

        bool write_solution = false;
        bool write_error = false;

        bool io_write_solution = false;
        bool io_write_error = false;

        std::vector<int> index_level;
        std::map<T_Key, int> map;



        //--------------------------------------------------------------------------------------------------------------

        BraidHeatCheck() : ::ug::xbraid::IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {
            index_level = std::vector<int>();
            err_u = BiotErrorData<TDomain, TAlgebra>();
            err_sol = BiotErrorData<TDomain, TAlgebra>();
            err_udiffsol = BiotErrorData<TDomain, TAlgebra>();
        };

        ~BraidHeatCheck() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_log(SP_ParallelLogger log) {
            this->m_log = log;
        }

        void set_solution_name(SP_VTKOutput vtk, const char *fname) {
            this->m_out_solution = make_sp(new T_VTK_ProcessObserver (vtk, fname));
            this->m_ioout_solution = make_sp(new T_IO_ProcessObserver (fname));

        }

        void set_diff_name(SP_VTKOutput vtk, const char *fname) {
            this->m_out_diff = make_sp(new T_VTK_ProcessObserver(vtk, fname));
            this->m_ioout_diff = make_sp(new T_IO_ProcessObserver(fname));
        }

        void set_num_ref(int ref) {
            this->num_ref = ref;
        }

        void set_max_index(int precomputed, int problem) {
            this->max_index_precomputed = precomputed;
            this->max_index = problem;
            index_level.resize(1);
            index_level[0] = this->max_index;
        }

        void set_c_factor(int level, int factor) {
            index_level.resize(level + 2);
            index_level[level + 1] = index_level[level] / factor;
            this->m_log->o << "level: " << index_level[level] << "\t" << factor << "\t" << index_level[level + 1]
                           << std::endl;
        }

        void set_vtk_write_mode(bool solution, bool error) {
            this->write_solution = solution;
            this->write_error = error;
        }

        void set_io_write_mode(bool solution, bool error) {
            this->io_write_solution = solution;
            this->io_write_error = error;
        }

        bool lua_write(SP_GridFunction u, int index, double time, double dt) {
            return this->step_process(u, index, time,0);
        }

        bool print(std::string str,SP_GridFunction u, int index, double time) {
            return this->step_process(u, index, time,0.0);
        }

        bool step_process(SP_GridFunction u, int index, double time, double dt) override {

                SP_GridFunction sol = u->clone_without_values();
                SP_GridFunction udiffsol = u->clone();

                // write vtk output
                if (this->write_solution) {
                    m_out_solution->step_process(u, index, time,0);
                }

                if (this->io_write_solution) {
                    m_ioout_solution->step_process(u, index, time,0);
                }


                // load gridfunction file (ref solution)

                Interpolate("AnalyticSolution", sol, "t",time);

                // substract
                ::ug::VecSubtract(*udiffsol.get(), *udiffsol.get(), *sol.get());

                // write vtk error
                if (this->write_error) {
                    m_out_diff->step_process(udiffsol, index, time,dt);
                }

                if (this->io_write_error) {
                    m_ioout_diff->step_process(udiffsol, index, time,dt);
                }

                // write norms
                m_log->o << ">> norms idx=" << index
                     << " t=" << time
                     << " iter=" << 0
                     << " level=" << 0
                     << " c=" << 0
                     << " done=" << 1
                     << std::endl;

                m_log->o << std::setw(10) << ">> norm"
                         << std::setw(20) << "solution"
                         << std::setw(20) << "error"
                         << std::setw(20) << "relative"
                         << std::endl;

            double nrm_u = L2Norm(*u.get(), "t", temporder);
            double nrm_sol =  L2Norm(*sol.get(), "t", temporder);
            double nrm_diff = L2Norm(*udiffsol.get(), "t", temporder);
            if (nrm_sol < 0.003) {
                nrm_sol = 0.003;
            }

            double nrm_rel = nrm_diff / nrm_sol;
            m_log->o << std::setw(10) << ">> l2(t)"
                     << std::setw(20) << nrm_u
                     << std::setw(20) << nrm_diff
                     << std::setw(20) << nrm_rel
                     << std::endl;

            return false; // no error
        };

        int temporder = 2 ;

        bool step_process(SP_GridFunction u, int index, double time, double dt, int iteration, int level) override {

            int count = 0;
            auto tuple = std::make_tuple(index, iteration, level);
            auto it = map.find(tuple);
            if (it != map.end()) {
                count = it->second;
                count += 1;
                map[tuple] = count;
            } else {
                count = 0;
                map.emplace(tuple, 0);
            }



                SP_GridFunction sol = u->clone_without_values();
                SP_GridFunction udiffsol = u->clone();

                // write vtk output
                if (this->write_solution) {
                    m_out_solution->step_process(u, index, time,0, iteration, level);
                }
                if (this->io_write_solution) {
                    m_ioout_solution->step_process(u, index, time, 0,iteration, level);
                }



                // load gridfunction file (ref solution)
                Interpolate("AnalyticSolution", sol, "t",time);





                // substract
                ::ug::VecSubtract(*udiffsol.get(),*udiffsol.get(), *sol.get());

                // write vtk error
                if (this->write_solution) {
                    m_out_diff->step_process(udiffsol, index, time, dt,iteration, level);
                }
                if (this->io_write_solution) {
                    m_ioout_diff->step_process(udiffsol, index, time, dt,iteration, level);
                }

                m_log->o << ">> norms idx=" << index
                     << " t=" << time
                     << " iter=" << iteration
                     << " level=" << level
                     << " c=" << count
                     << " done=" << false
                     << std::endl;

                m_log->o << std::setw(10) << ">> norm"
                         << std::setw(20) << "solution"
                         << std::setw(20) << "error"
                         << std::setw(20) << "relative"
                         << std::endl;

                double nrm_u = L2Norm(*u.get(), "t", temporder);
                double nrm_sol =  L2Norm(*sol.get(), "t", temporder);
                double nrm_diff = L2Norm(*udiffsol.get(), "t", temporder);
                if (nrm_sol < 0.003) {
                    nrm_sol = 0.003;
                }

            double nrm_rel = nrm_diff / nrm_sol;
                m_log->o << std::setw(10) << ">> l2(t)"
                         << std::setw(20) << nrm_u
                         << std::setw(20) << nrm_diff
                         << std::setw(20) << nrm_rel
                         << std::endl;



            return false; // no error
        };



        bool lua_compare(SP_GridFunction u,SP_GridFunction v, int index, double time, int iteration, int level) {
            return false;
        };
        //--------------------------------------------------------------------------------------------------------------
    };
}}}

#endif