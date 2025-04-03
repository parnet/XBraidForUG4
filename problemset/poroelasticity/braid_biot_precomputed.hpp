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

        BraidBiotCheckPrecomputed() : IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {
            index_level_ = std::vector<int>();
            err_u_ = BiotErrorData<TDomain, TAlgebra>();
            err_sol_ = BiotErrorData<TDomain, TAlgebra>();
            err_udiffsol_ = BiotErrorData<TDomain, TAlgebra>();
        };

        ~BraidBiotCheckPrecomputed() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_log(SP_ParallelLogger log) {
            this->log_ = log;
        }

        void set_base_path(std::string path) {
            this->base_path_ = path;
        }

        void set_solution_name(SP_VTKOutput vtk, const char *fname) {
            this->out_solution_ = make_sp(new T_VTK_ProcessObserver (vtk, fname));
            this->ioout_solution_ = make_sp(new T_IO_ProcessObserver (fname));

        }

        void set_diff_name(SP_VTKOutput vtk, const char *fname) {
            this->out_diff_ = make_sp(new T_VTK_ProcessObserver(vtk, fname));
            this->ioout_diff_ = make_sp(new T_IO_ProcessObserver(fname));
        }

        void compare_norms(int index, double time, int iteration, int level, int c, bool done) {
            log_->o << ">> norms idx=" << index
                     << " t=" << time
                     << " iter=" << iteration
                     << " level=" << level
                     << " c=" << c
                     << " done=" << done
                     << std::endl;

            log_->o << std::setw(10) << ">> norm"
                     << std::setw(20) << "solution"
                     << std::setw(20) << "error"
                     << std::setw(20) << "relative"
                     << std::endl;

            log_->o << std::setw(10) << ">> l2(p)"
                     << std::setw(20) << err_u_.l2_norm_p_
                     << std::setw(20) << err_udiffsol_.l2_norm_p_
                     << std::setw(20) << (err_udiffsol_.l2_norm_p_ / err_sol_.l2_norm_p_)
                     << std::endl;

            log_->o << std::setw(10) << ">> l2(ux)"
                     << std::setw(20) << err_u_.l2_norm_ux_
                     << std::setw(20) << err_udiffsol_.l2_norm_ux_
                     << std::setw(20) << (err_udiffsol_.l2_norm_ux_ / err_sol_.l2_norm_ux_)
                     << std::endl;

            log_->o << std::setw(10) << ">> l2(uy)"
                     << std::setw(20) << err_u_.l2_norm_uy_
                     << std::setw(20) << err_udiffsol_.l2_norm_uy_
                     << std::setw(20) << (err_udiffsol_.l2_norm_uy_ / err_sol_.l2_norm_uy_)
                     << std::endl;


            log_->o << std::setw(10) << ">> h1(ux)"
                     << std::setw(20) << err_u_.h1_norm_ux_
                     << std::setw(20) << err_udiffsol_.h1_norm_ux_
                     << std::setw(20) << (err_udiffsol_.h1_norm_ux_ / err_sol_.h1_norm_ux_)
                     << std::endl;

            log_->o << std::setw(10) << ">> h1(uy)"
                     << std::setw(20) << err_u_.h1_norm_uy_
                     << std::setw(20) << err_udiffsol_.h1_norm_uy_
                     << std::setw(20) << (err_udiffsol_.h1_norm_uy_ / err_sol_.h1_norm_uy_)
                     << std::endl << std::endl;
        }

        void set_num_ref(int ref) {
            this->num_ref_ = ref;
        }

        void set_max_index(int precomputed, int problem) {
            this->max_index_precomputed_ = precomputed;
            this->max_index_ = problem;
            index_level_.resize(1);
            index_level_[0] = this->max_index_;
        }

        void set_c_factor(int level, int factor) {
            index_level_.resize(level + 2);
            index_level_[level + 1] = index_level_[level] / factor;
            this->log_->o << "level: " << index_level_[level] << "\t" << factor << "\t" << index_level_[level + 1]
                           << std::endl;
        }

        void set_vtk_write_mode(bool solution, bool error) {
            this->write_solution_ = solution;
            this->write_error_ = error;
        }

        void set_io_write_mode(bool solution, bool error) {
            this->io_write_solution_ = solution;
            this->io_write_error_ = error;
        }

        bool lua_write(SP_GridFunction u, int index, double time, double dt) {
            return this->step_process(u, index, time,0);
        }

        bool print(std::string str,SP_GridFunction u, int index, double time) {
            return this->step_process(u, index, time,0.0);
        }

        bool step_process(SP_GridFunction u, int index, double time, double dt) override {
            int zidx = (index * this->max_index_precomputed_) / index_level_[0];
            int rem = (index * this->max_index_precomputed_) % index_level_[0];
            if (rem == 0 && index != 0) {
                SP_GridFunction sol = u->clone_without_values();
                SP_GridFunction udiffsol = u->clone();

                // write vtk output
                if (this->write_solution_) {
                    out_solution_->step_process(u, index, time,0);
                }

                if (this->io_write_solution_) {
                    ioout_solution_->step_process(u, index, time,0);
                }


                // load gridfunction file (ref solution)
                xbraid::PIOGridFunction <TDomain, TAlgebra> io = xbraid::PIOGridFunction<TDomain, TAlgebra>();
                std::stringstream ss_ref;
                int local_num_ref = u->grid_level().level();
                if(local_num_ref == GridLevel::TOP) {
                    local_num_ref = this->num_ref_;
                }
                ss_ref << this->base_path_ << "/num_ref_" << local_num_ref << "/BarryMercer_2D_NumRef"<< local_num_ref<<"_nX1_"<< zidx;
                io.read(sol, ss_ref.str().c_str());

                // substract
                ::ug::VecSubtract(*udiffsol.get(), *udiffsol.get(), *sol.get());

                // write vtk error
                if (this->write_error_) {
                    out_diff_->step_process(udiffsol, index, time,dt);
                }

                if (this->io_write_error_) {
                    ioout_diff_->step_process(udiffsol, index, time,dt);
                }

                // compute norms
                err_u_.compute(u->clone());
                err_sol_.compute(sol->clone());
                err_udiffsol_.compute(udiffsol->clone());

                // write norms
                compare_norms(index, time, 0, 0, 0, true);
            }
            return false; // no error
        };

        bool step_process(SP_GridFunction u, int index, double time, double dt, int iteration, int level) override {

            int count = 0;
            auto tuple = std::make_tuple(index, iteration, level);
            auto it = map_.find(tuple);
            if (it != map_.end()) {
                count = it->second;
                count += 1;
                map_[tuple] = count;
            } else {
                count = 0;
                map_.emplace(tuple, 0);
            }


            int zidx = (index * this->max_index_precomputed_) / index_level_[level];
            int rem = (index * this->max_index_precomputed_) % index_level_[level];

            if (rem == 0 && index != 0) {
                SP_GridFunction sol = u->clone_without_values();
                SP_GridFunction udiffsol = u->clone();

                // write vtk output
                if (this->write_solution_) {
                    out_solution_->step_process(u, index, time,0, iteration, level);
                }
                if (this->io_write_solution_) {
                    ioout_solution_->step_process(u, index, time, 0,iteration, level);
                }

                // load gridfunction file (ref solution)
                PIOGridFunction <TDomain, TAlgebra> io = xbraid::PIOGridFunction<TDomain, TAlgebra>();
                std::stringstream ss_ref;
                int local_num_ref = u->grid_level().level();
                if(local_num_ref == GridLevel::TOP) {
                    local_num_ref = this->num_ref_;
                }
                ss_ref << this->base_path_ << "/num_ref_" << local_num_ref << "/BarryMercer_2D_NumRef"<< local_num_ref<<"_nX1_"<< zidx;
                std::cout << index << "\t";
                std::cout << ss_ref.str().c_str() << std::endl;
                io.read(sol, ss_ref.str().c_str());
                // substract
                ::ug::VecSubtract(*udiffsol.get(),*udiffsol.get(), *sol.get());

                // write vtk error
                if (this->write_solution_) {
                    out_diff_->step_process(udiffsol, index, time, dt,iteration, level);
                }
                if (this->io_write_solution_) {
                    ioout_diff_->step_process(udiffsol, index, time, dt,iteration, level);
                }

                // compute norms
                err_u_.compute(u->clone());
                err_sol_.compute(sol->clone());
                err_udiffsol_.compute(udiffsol->clone());

                // write norms
                compare_norms(index, time, iteration, level, count, false);
            }

            return false; // no error
        };



        bool lua_compare(SP_GridFunction u,SP_GridFunction v, int index, double time, int iteration, int level) {

            int count = 0;
            auto tuple = std::make_tuple(index, iteration, level);
            auto it = map_.find(tuple);
            if (it != map_.end()) {
                count = it->second;
                count += 1;
                map_[tuple] = count;
            } else {
                count = 0;
                map_.emplace(tuple, 0);
            }


            SP_GridFunction udiffsol = u->clone();

            // write vtk output
            if (this->write_solution_) {
                out_solution_->step_process(u, index, time, 0,iteration, level);
            }
            if (this->io_write_solution_) {
                ioout_solution_->step_process(u, index, time, 0,iteration, level);
            }

            // load gridfunction file (ref solution)
            ::ug::VecSubtract(*udiffsol.get(), *udiffsol.get(), *v.get());

            // write vtk error
            if (this->write_solution_) {
                out_diff_->step_process(udiffsol, index, time,0, iteration, level);
            }
            if (this->io_write_solution_) {
                ioout_diff_->step_process(udiffsol, index, time,0, iteration, level);
            }

            // compute norms
            err_u_.compute(u->clone());
            err_sol_.compute(v->clone());
            err_udiffsol_.compute(udiffsol->clone());

            // write norms
            compare_norms(index, time, iteration, level, count, false);
            return false; // no error
        };
        //--------------------------------------------------------------------------------------------------------------

        BiotErrorData<TDomain, TAlgebra> err_u_;
        BiotErrorData<TDomain, TAlgebra> err_sol_;
        BiotErrorData<TDomain, TAlgebra> err_udiffsol_;

        SP_VTK_ProcessObserver out_solution_;
        SP_VTK_ProcessObserver out_diff_;

        SP_IO_ProcessOobserver ioout_solution_;
        SP_IO_ProcessOobserver ioout_diff_;

        SP_ParallelLogger log_;

        int num_ref_ = 3;
        int max_index_ = 512;
        int max_index_precomputed_ = 512;

        bool write_solution_ = false;
        bool write_error_ = false;

        bool io_write_solution_ = false;
        bool io_write_error_ = false;

        std::string base_path_ = "/storage/data/poroelasticity/analyticsolution";
        std::vector<int> index_level_;
        std::map<T_Key, int> map_;

        //--------------------------------------------------------------------------------------------------------------
    };
}}}
#endif