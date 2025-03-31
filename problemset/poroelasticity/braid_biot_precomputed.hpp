#ifndef UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRAID_BIOT_PRECOMPUTED_HPP
#define UGPLUGIN_XBRAIDFORUG4_POROELASTICITY_BRAID_BIOT_PRECOMPUTED_HPP

#include "Poroelasticity/biot_tools.h"

#include "interface/observer_xbraid.hpp"
#include "observer/io_process_observer.hpp"
#include "observer/vtk_process_observer.hpp"
#include "util/parallel_logger.hpp"
#include "util/parallel_io_gridfunction.hpp"

#include "braid_biot_estimator.hpp"
#include "biot_error_data.hpp"
//2025-03 #include <cstdint>

namespace ug { namespace xbraid{ namespace poro {

    template<typename TDomain, typename TAlgebra>
    class BraidBiotCheckPrecomputed : public IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
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

        std::string base_path = "/storage/data/poroelasticity/analyticsolution";
        std::vector<int> index_level;
        std::map<T_Key, int> map;



        //--------------------------------------------------------------------------------------------------------------

        BraidBiotCheckPrecomputed() : IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {
            index_level = std::vector<int>();
            err_u = BiotErrorData<TDomain, TAlgebra>();
            err_sol = BiotErrorData<TDomain, TAlgebra>();
            err_udiffsol = BiotErrorData<TDomain, TAlgebra>();
        };

        ~BraidBiotCheckPrecomputed() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_log(SP_ParallelLogger log) {
            this->m_log = log;
        }

        void set_base_path(std::string path) {
            this->base_path = path;
        }

        void set_solution_name(SP_VTKOutput vtk, const char *fname) {
            this->m_out_solution = make_sp(new T_VTK_ProcessObserver (vtk, fname));
            this->m_ioout_solution = make_sp(new T_IO_ProcessObserver (fname));

        }

        void set_diff_name(SP_VTKOutput vtk, const char *fname) {
            this->m_out_diff = make_sp(new T_VTK_ProcessObserver(vtk, fname));
            this->m_ioout_diff = make_sp(new T_IO_ProcessObserver(fname));
        }

        void compare_norms(int index, double time, int iteration, int level, int c, bool done) {
            m_log->o << ">> norms idx=" << index
                     << " t=" << time
                     << " iter=" << iteration
                     << " level=" << level
                     << " c=" << c
                     << " done=" << done
                     << std::endl;

            m_log->o << std::setw(10) << ">> norm"
                     << std::setw(20) << "solution"
                     << std::setw(20) << "error"
                     << std::setw(20) << "relative"
                     << std::endl;

            m_log->o << std::setw(10) << ">> l2(p)"
                     << std::setw(20) << err_u.l2_norm_p
                     << std::setw(20) << err_udiffsol.l2_norm_p
                     << std::setw(20) << (err_udiffsol.l2_norm_p / err_sol.l2_norm_p)
                     << std::endl;

            m_log->o << std::setw(10) << ">> l2(ux)"
                     << std::setw(20) << err_u.l2_norm_ux
                     << std::setw(20) << err_udiffsol.l2_norm_ux
                     << std::setw(20) << (err_udiffsol.l2_norm_ux / err_sol.l2_norm_ux)
                     << std::endl;

            m_log->o << std::setw(10) << ">> l2(uy)"
                     << std::setw(20) << err_u.l2_norm_uy
                     << std::setw(20) << err_udiffsol.l2_norm_uy
                     << std::setw(20) << (err_udiffsol.l2_norm_uy / err_sol.l2_norm_uy)
                     << std::endl;


            m_log->o << std::setw(10) << ">> h1(ux)"
                     << std::setw(20) << err_u.h1_norm_ux
                     << std::setw(20) << err_udiffsol.h1_norm_ux
                     << std::setw(20) << (err_udiffsol.h1_norm_ux / err_sol.h1_norm_ux)
                     << std::endl;

            m_log->o << std::setw(10) << ">> h1(uy)"
                     << std::setw(20) << err_u.h1_norm_uy
                     << std::setw(20) << err_udiffsol.h1_norm_uy
                     << std::setw(20) << (err_udiffsol.h1_norm_uy / err_sol.h1_norm_uy)
                     << std::endl << std::endl;
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
            int zidx = (index * this->max_index_precomputed) / index_level[0];
            int rem = (index * this->max_index_precomputed) % index_level[0];
            if (rem == 0 && index != 0) {
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
                xbraid::PIOGridFunction <TDomain, TAlgebra> io = xbraid::PIOGridFunction<TDomain, TAlgebra>();
                std::stringstream ss_ref;
                int local_num_ref = u->grid_level().level();
                if(local_num_ref == GridLevel::TOP) {
                    local_num_ref = this->num_ref;
                }
                ss_ref << this->base_path << "/num_ref_" << local_num_ref << "/BarryMercer_2D_NumRef"<< local_num_ref<<"_nX1_"<< zidx;
                io.read(sol, ss_ref.str().c_str());

                // substract
                ::ug::VecSubtract(*udiffsol.get(), *udiffsol.get(), *sol.get());

                // write vtk error
                if (this->write_error) {
                    m_out_diff->step_process(udiffsol, index, time,dt);
                }

                if (this->io_write_error) {
                    m_ioout_diff->step_process(udiffsol, index, time,dt);
                }

                // compute norms
                err_u.compute(u->clone());
                err_sol.compute(sol->clone());
                err_udiffsol.compute(udiffsol->clone());

                // write norms
                compare_norms(index, time, 0, 0, 0, true);
            }
            return false; // no error
        };

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


            int zidx = (index * this->max_index_precomputed) / index_level[level];
            int rem = (index * this->max_index_precomputed) % index_level[level];

            if (rem == 0 && index != 0) {
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
                PIOGridFunction <TDomain, TAlgebra> io = xbraid::PIOGridFunction<TDomain, TAlgebra>();
                std::stringstream ss_ref;
                int local_num_ref = u->grid_level().level();
                if(local_num_ref == GridLevel::TOP) {
                    local_num_ref = this->num_ref;
                }
                ss_ref << this->base_path << "/num_ref_" << local_num_ref << "/BarryMercer_2D_NumRef"<< local_num_ref<<"_nX1_"<< zidx;
                std::cout << index << "\t";
                std::cout << ss_ref.str().c_str() << std::endl;
                io.read(sol, ss_ref.str().c_str());
                // substract
                ::ug::VecSubtract(*udiffsol.get(),*udiffsol.get(), *sol.get());

                // write vtk error
                if (this->write_solution) {
                    m_out_diff->step_process(udiffsol, index, time, dt,iteration, level);
                }
                if (this->io_write_solution) {
                    m_ioout_diff->step_process(udiffsol, index, time, dt,iteration, level);
                }

                // compute norms
                err_u.compute(u->clone());
                err_sol.compute(sol->clone());
                err_udiffsol.compute(udiffsol->clone());

                // write norms
                compare_norms(index, time, iteration, level, count, false);
            }

            return false; // no error
        };



        bool lua_compare(SP_GridFunction u,SP_GridFunction v, int index, double time, int iteration, int level) {

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


            SP_GridFunction udiffsol = u->clone();

            // write vtk output
            if (this->write_solution) {
                m_out_solution->step_process(u, index, time, 0,iteration, level);
            }
            if (this->io_write_solution) {
                m_ioout_solution->step_process(u, index, time, 0,iteration, level);
            }

            // load gridfunction file (ref solution)
            ::ug::VecSubtract(*udiffsol.get(), *udiffsol.get(), *v.get());

            // write vtk error
            if (this->write_solution) {
                m_out_diff->step_process(udiffsol, index, time,0, iteration, level);
            }
            if (this->io_write_solution) {
                m_ioout_diff->step_process(udiffsol, index, time,0, iteration, level);
            }

            // compute norms
            err_u.compute(u->clone());
            err_sol.compute(v->clone());
            err_udiffsol.compute(udiffsol->clone());

            // write norms
            compare_norms(index, time, iteration, level, count, false);
            return false; // no error
        };
        //--------------------------------------------------------------------------------------------------------------
    };
}}}
#endif