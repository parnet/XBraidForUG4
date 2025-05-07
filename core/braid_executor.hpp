#ifndef UUGPLUGIN_XBRAIDFORUG4_CORE_BRAID_EXECUTOR_H
#define UGPLUGIN_XBRAIDFORUG4_CORE_BRAID_EXECUTOR_H

#include "bindings/lua/lua_user_data.h"
#include "driver/basic_driver.hpp"
#include "interface/observer_xbraid.hpp"
#include "space_time_communicator.hpp"
#include "braid_settings.hpp"
#include "initializer/start_value_initializer.hpp"

#include "config/pragma.hpp"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidExecutor {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_BraidBaseDriver = BraidGridFunctionBase<TDomain, TAlgebra> ;
        using SP_BraidBaseDriver = SmartPtr<T_BraidBaseDriver> ;

        using SP_SpaceTimeCommunicator = SmartPtr<SpaceTimeCommunicator> ;

        using T_ITimeIntegratorObserver = ITimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_T_ITimeIntegratorObserver = SmartPtr<T_ITimeIntegratorObserver> ;

        using T_IXBraidTimeIntegratorObserver = IXBraidTimeIntegratorObserver<TDomain, TAlgebra> ;
        using SP_XBTIObserver = SmartPtr<T_IXBraidTimeIntegratorObserver> ;

        using T_BraidInitializer = BraidInitializer<TDomain, TAlgebra> ;
        using SP_BraidInitializer = SmartPtr<T_BraidInitializer> ;

        using T_ParallelLogger = ParallelLogger ;
        using SP_ParallelLogger = SmartPtr<T_ParallelLogger> ;

        using TSpatialNorm = BraidSpatialNorm<TDomain, TAlgebra> ;
        using SPSpatialNorm = SmartPtr<TSpatialNorm> ;



        //--------------------------------------------------------------------------------------------------------------

        BraidExecutor(SP_SpaceTimeCommunicator comm, SP_BraidBaseDriver driver) {
            this->comm_ = comm;
            this->set_driver(driver);
        }

        virtual ~BraidExecutor() {
            delete braid_core_;
        }

        //--------------------------------------------------------------------------------------------------------------

        void set_residual(bool residual) {
            if (residual ) { // && this->m_driver->can_residual_method
                this->braid_settings_.residual_ = residual;
                this->braid_core_->SetResidual();
            }
        }

        void set_n_relax(int level, int number) {
            if (level == -1) {
                this->braid_settings_.n_relax_default_ = number;
            }
            else {
                this->braid_settings_.n_relax_[level] = number;
                if (level > this->braid_settings_.eval_level_) {
                    this->braid_settings_.eval_level_ = level;
                }
            }
            this->braid_core_->SetNRelax(level, number);
        }

        void set_c_factor(int level, int factor) {
            if (level == -1) {
                this->braid_settings_.c_factor_default_ = factor;
            }
            else {
                this->braid_settings_.c_factor_[level] = factor;
                if (level > this->braid_settings_.eval_level_) {
                    this->braid_settings_.eval_level_ = level;
                }
            }
            this->braid_core_->SetCFactor(level, factor);
        }

        void set_max_levels(int maxLevel) {
            this->braid_settings_.max_level_ = maxLevel;
            this->braid_core_->SetMaxLevels(maxLevel);
            this->driver_->set_max_levels(maxLevel);
        }

        void set_skip_downcycle_work(bool skip) {
            this->braid_settings_.skip_ = skip;
            this->braid_core_->SetSkip(skip);
        }

        void set_min_coarse(int minCoarse) {
            this->braid_settings_.min_coarse_ = minCoarse;
            this->braid_core_->SetMinCoarse(minCoarse);
        }

        void set_max_iterations(int max_iter) {
            this->braid_settings_.max_iter_ = max_iter;
            this->braid_core_->SetMaxIter(max_iter);
        }

        void set_absolute_tol(double tol) {
            this->braid_settings_.abs_tol_ = tol;
            this->braid_core_->SetAbsTol(tol);
        }

        void set_relative_tol(double tol) {
            this->braid_settings_.rel_tol_ = tol;
            this->braid_core_->SetRelTol(tol);
        }

        void set_temporal_norm(int nrm) {
            this->braid_settings_.temp_norm_ = nrm;
            this->braid_core_->SetTemporalNorm(nrm);
        }

        void set_sequential(bool sequential) {
            this->braid_settings_.sequential_ = sequential;
            if (sequential) {
                this->braid_core_->SetSeqSoln(1);
            } else {
                this->braid_core_->SetSeqSoln(0);
            }
        }

        void set_store_values(int level) {
            this->braid_settings_.store_values_ = level;
            this->braid_core_->SetStorage(level);
        }

        void set_spatial_coarsen_and_refine(bool cnr) {
            if (cnr) {
                this->braid_settings_.coarsen_and_refine_ = cnr;
                this->braid_core_->SetSpatialCoarsenAndRefine();
            }
        }

        void set_refine(bool ref) {
            this->braid_settings_.refine_ = ref;
            if (ref) {
                this->braid_core_->SetRefine(1);
            } else {
                this->braid_core_->SetRefine(0);
            }
        }

        void set_max_refinements(int number) {
            this->braid_settings_.max_refinements_ = number;
            this->braid_core_->SetMaxRefinements(number);
        }

        void set_access_level(int level) {
            this->braid_settings_.access_level_ = level;
            this->braid_core_->SetAccessLevel(level);
        }

        void set_print_level(int level) {
            this->braid_settings_.print_level_ = level;
            this->braid_core_->SetPrintLevel(level);
        }

        void set_default_print_file() {
            const char* file = "braid_runtime.out";
            this->set_print_file(file);
        }

        void set_print_file(const char* file) {
            this->braid_settings_.print_file_ = file;
            this->braid_core_->SetPrintFile(file);
        }

        void set_cycle_fmg() {
            this->braid_settings_.cycle_fmg_ = true;
            this->braid_core_->SetFMG();
        }

        void set_cycle_nfmg(int vu) {
            this->braid_settings_.cycle_nfmg_ = vu;
            this->braid_core_->SetNFMG(vu);
        }

        void set_cycle_nfmgv(int mu) {
            this->braid_settings_.cycle_nfmgv_ = mu;
            this->braid_core_->SetNFMGVcyc(mu);
        }

        void set_cycle_type(const char* ctyp) {
            std::cout << "debBraidExecutor::setCycleType(args)" << std::endl << std::flush;
            if (strcmp(ctyp, "V FCF") == 0) {
                this->braid_core_->SetNRelax(-1, 1);
                this->braid_core_->SetNRelax(0, 1);
            } else if (strcmp(ctyp, "V F") == 0) {
                this->braid_core_->SetNRelax(-1, 0);
                this->braid_core_->SetNRelax(0, 0);
            } else if (strcmp(ctyp, "V F-FCF") == 0) {
                this->braid_core_->SetNRelax(-1, 1);
                this->braid_core_->SetNRelax(0, 0);
            } else if (strcmp(ctyp, "F FCF") == 0) {
                this->braid_core_->SetNRelax(-1, 1);
                this->braid_core_->SetNRelax(0, 1);
                this->braid_settings_.cycle_fmg_ = true;
                this->braid_core_->SetFMG();
            } else if (strcmp(ctyp, "F F") == 0) {
                this->braid_core_->SetNRelax(-1, 0);
                this->braid_core_->SetNRelax(0, 0);
                this->braid_settings_.cycle_fmg_ = true;
                this->braid_core_->SetFMG();
            } else if (strcmp(ctyp, "F F-FCF") == 0) {
                this->braid_core_->SetNRelax(-1, 1);
                this->braid_core_->SetNRelax(0, 0);
                this->braid_settings_.cycle_fmg_ = true;
                this->braid_core_->SetFMG();
            } else {
                std::cout << "invalid cycle type" << std::endl;
            }
        }

        void set_sync() {
            this->braid_settings_.sync_ = true;
            this->braid_core_->SetSync();
        }

        void set_increase_max_levels() {
            this->braid_settings_.increase_max_level_ = true;
            this->braid_core_->SetIncrMaxLevels();
        }

        void set_relax_only_cg(bool setting) {
            this->braid_settings_.relax_only_cg_ = setting;
            if (setting) {
                this->braid_core_->SetRelaxOnlyCG(1);
            } else {
                this->braid_core_->SetRelaxOnlyCG(0);
            }
        }

        void set_agg_c_factor(int cfactor0) {
            this->braid_core_->SetAggCFactor(cfactor0);
        }

        void set_periodic(int periodic) {
            this->braid_settings_.periodic_ = periodic;
            this->braid_core_->SetPeriodic(periodic);
        }

        void set_final_fc_relax() {
            this->braid_settings_.final_fc_relax_ = true;
            this->braid_core_->SetFinalFCRelax();
        }

        void set_reverted_ranks(int ranks) {
            this->braid_settings_.reverted_ranks_ = ranks;
            this->braid_core_->SetRevertedRanks(ranks);
        }

        void set_richardson_estimation(bool use_richardson, bool use_extrapolation, int local_order) {
            this->braid_settings_.est_error_ = use_richardson;
            this->braid_settings_.richardson_ = use_extrapolation;
            this->braid_settings_.local_order_ = local_order;
            this->braid_core_->SetRichardsonEstimation(use_richardson, use_extrapolation, local_order);
        }

        void set_file_io_level(int level) {
            this->braid_settings_.file_io_level_ = level;
            this->braid_core_->SetFileIOLevel(level);
        }

        void set_c_relax_weight(int level, double weight) {
            if (level == -1) {
                this->braid_settings_.c_relax_weight_default_ = weight;
            } else {
                this->braid_settings_.c_relax_weight_[level] = weight;
                if (level > this->braid_settings_.eval_level_) {
                    this->braid_settings_.eval_level_ = level;
                }
            }
            this->braid_core_->SetCRelaxWt(level, weight);
        }

        void set_t_points_cutoff(int cutoff) {
            this->braid_settings_.t_points_cutoff_ = cutoff;
            this->braid_core_->SetTPointsCutoff(cutoff);
        }

        void set_full_residual_norm() {
            this->braid_settings_.full_r_norm_ = true;
            this->braid_core_->SetFullRNormRes(_BraidAppResidual);
        }

        void set_time_grid() {}

        int get_num_iteration() {
            int iter = 0;
            this->braid_core_->GetNumIter(&iter);
            return iter;
        }

        void get_c_factor() {}

        void get_residual_norms() {}

        int get_num_level() {
            int number_of_level = 0;
            this->braid_core_->GetNLevels(&number_of_level);
            return number_of_level;
        }

        int get_warm_restart() const {
            return this->braid_core_->GetWarmRestart();
        }

        int get_distribution_lower() const {
            int lower;
            int upper;
            this->braid_core_->GetDistribution(&lower, &upper);
            return lower;
        }

        int get_distribution_upper() const {
            int lower;
            int upper;
            this->braid_core_->GetDistribution(&lower, &upper);
            return upper;
        }

        void get_distribution(int& lower, int& upper) const {
            this->braid_core_->GetDistribution(&lower, &upper);
        }

        int get_id() const {
            int id = 0;
            this->braid_core_->GetMyID(&id);
            return id;
        }

        braid_Core get_real_braid_core() const {
            return this->braid_core_->GetCore();
        }

        BraidCore* get_braid_core() const {
            return this->braid_core_;
        }

        void reset_core() {
            std::cout << "core_ is " << braid_core_ << std::endl;
            if (braid_core_ != nullptr) {
                delete braid_core_;
            }
            braid_core_ = new BraidCore(this->comm_->GLOBAL, this->driver_.get());
        }

        void set_driver(SP_BraidBaseDriver p_app) {
            this->driver_ = p_app;
            this->driver_->comm_ = this->comm_;
            this->driver_->comm_t = this->comm_->TEMPORAL;  // sets the temporal communicator for the library
            this->reset_core();
        }

        SP_BraidBaseDriver get_driver() {
            return this->driver_;
        }

        void print_settings() {
            // todo move to settings & make parseable
            this->log_->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
            this->log_->o << "Module Name:                                             Braid Executor" << std::endl;
            this->log_->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----" << std::endl;
            this->log_->o << "Max Level: " << this->braid_settings_.max_level_ << std::endl;
            this->log_->o << "Increase max level: " << this->braid_settings_.increase_max_level_ << std::endl;
            this->log_->o << "Skip Down Cylce Work: " << this->braid_settings_.skip_ << std::endl;
            this->log_->o << "Absolute Tol: " << this->braid_settings_.abs_tol_ << std::endl;
            this->log_->o << "Relative Tol: " << this->braid_settings_.rel_tol_ << std::endl;
            this->log_->o << "Temporal Norm: " << this->braid_settings_.temp_norm_ << std::endl;
            this->log_->o << "Max Iteration: " << this->braid_settings_.max_iter_ << std::endl;
            this->log_->o << "Sync: " << this->braid_settings_.sync_ << std::endl;
            this->log_->o << "Residual: " << this->braid_settings_.residual_ << std::endl;
            this->log_->o << "Cycle FMG: " << this->braid_settings_.cycle_fmg_ << std::endl;
            this->log_->o << "Cycle NFMG: " << this->braid_settings_.cycle_nfmg_ << std::endl;
            this->log_->o << "Cycle NFMGV: " << this->braid_settings_.cycle_nfmgv_ << std::endl;
            this->log_->o << "Store Values: " << this->braid_settings_.store_values_ << std::endl;
            this->log_->o << "Print File: " << this->braid_settings_.print_file_ << std::endl;
            this->log_->o << "Print Level: " << this->braid_settings_.print_level_ << std::endl;
            this->log_->o << "File IO Level: " << this->braid_settings_.file_io_level_ << std::endl;
            this->log_->o << "Access Level: " << this->braid_settings_.access_level_ << std::endl;
            this->log_->o << "Sequential: " << this->braid_settings_.sequential_ << std::endl;
            this->log_->o << "Coarsen & Refine: " << this->braid_settings_.coarsen_and_refine_ << std::endl;
            this->log_->o << "Min Coarse: " << this->braid_settings_.min_coarse_ << std::endl;
            this->log_->o << "Final FC Relax: " << this->braid_settings_.final_fc_relax_ << std::endl;
            this->log_->o << "Refine: " << this->braid_settings_.refine_ << std::endl;
            this->log_->o << "Max Refinements: " << this->braid_settings_.max_refinements_ << std::endl;
            this->log_->o << "Relax only CG: " << this->braid_settings_.relax_only_cg_ << std::endl;
            this->log_->o << "t-point cutoff: " << this->braid_settings_.t_points_cutoff_ << std::endl;
            this->log_->o << "Braid ID: " << this->get_id() << std::endl;
            this->log_->o << "Distribution: " << this->get_distribution_lower() << " , " << this->get_distribution_upper() << std::endl;
            this->log_->o << "Number of Level: " << this->get_num_level() << std::endl;
            this->log_->o << "Number of Iteration: " << this->get_num_iteration() << std::endl;
            this->log_->o << "Richardson-based Error Estimation: " << this->braid_settings_.est_error_ << std::endl;
            this->log_->o << "Richardson-based Extrapolation: " << this->braid_settings_.richardson_ << std::endl;
            this->log_->o << "Local order: " << this->braid_settings_.local_order_ << std::endl;
            this->log_->o << "Reverted Ranks: " << this->braid_settings_.reverted_ranks_ << std::endl;
            this->log_->o << "Periodic: " << this->braid_settings_.periodic_ << std::endl;
            this->log_->o << "Full residual norm: " << this->braid_settings_.full_r_norm_ << std::endl;
            this->log_->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----" << std::endl;
            this->log_->o << "default" << braid_settings_.n_relax_default_ << "    " << braid_settings_.c_relax_weight_default_ << "    " << braid_settings_.c_factor_default_ << "    " << std::endl;
            for (int i = 0; i <= this->braid_settings_.eval_level_; i++) {
                this->log_->o << i << "    ";
                if (this->braid_settings_.n_relax_.find(i) != this->braid_settings_.n_relax_.end()) {
                    this->log_->o << " " << this->braid_settings_.n_relax_[i] << " ";
                } else {
                    this->log_->o << "(" << braid_settings_.n_relax_default_ << ")";
                }
                this->log_->o << "    ";
                if (this->braid_settings_.c_relax_weight_.find(i) !=
                    this->braid_settings_.c_relax_weight_.end()) {
                    this->log_->o << " " << this->braid_settings_.c_relax_weight_[i] << " ";
                } else {
                    this->log_->o << "(" << braid_settings_.c_relax_weight_default_ << ")";
                }
                this->log_->o << "    ";
                if (this->braid_settings_.c_factor_.find(i) != this->braid_settings_.c_factor_.end()) {
                    this->log_->o << " " << this->braid_settings_.c_factor_[i] << " ";
                } else {
                    this->log_->o << "(" << braid_settings_.c_factor_default_ << ")";
                }
                this->log_->o << std::endl;
            }
            this->log_->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        }

        void print_summary() {
            // todo make parsable?
            this->log_->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----"
                << std::endl;
            this->log_->o << "Number of Level: " << this->get_num_level() << std::endl;
            this->log_->o << "Number of Iteration: " << this->get_num_iteration() << std::endl;
            int request = this->braid_settings_.max_iter_;
            auto* norms = new double[this->braid_settings_.max_iter_];
            this->braid_core_->GetRNorms(&request, norms);
            for (int i = 0; i < request; i++) {
                this->log_->o << i << ": " << norms[i] << std::endl;
            }
            this->log_->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        }

        void set_observer(SP_T_ITimeIntegratorObserver output) {
            this->driver_->set_observer(output);
        }

        void set_xbraid_observer(SP_XBTIObserver output) {
            this->driver_->set_xbraid_observer(output);
        }

        void set_initializer(SP_BraidInitializer initializer) {
            this->driver_->initializer_ = initializer;
        }

        void set_norm_provider(SPSpatialNorm norm) {
            this->driver_->norm_ = norm;
        }

        bool apply(SP_GridFunction u_stop, number t_stop, SP_GridFunction u0, number t0) {
            __debug(std::cout << "BraidExecutor::apply" << std::endl);

            this->driver_->set_start_time(t0);
            this->driver_->set_end_time(t_stop);
            this->driver_->set_start_vector(u0);

            if (this->driver_->initializer_ == SPNULL) {
                std::cout << "Create default Initializer " << std::endl;
                this->driver_->initializer_ = make_sp(new GridFunctionInitializer<TDomain,TAlgebra>());
            }

            if (this->log_ == SPNULL) {
                auto logger = make_sp(new T_ParallelLogger());
                logger->set_comm(this->comm_);
                logger->set_filename("std_output");
                logger->init();
                this->set_parallel_logger(logger);
            }

            this->driver_->initializer_->set_start_values(u0, t0);

            this->driver_->restimate_ = this->braid_settings_.est_error_;
            this->driver_->refine_time_ = this->braid_settings_.refine_;
            this->driver_->init();
            this->braid_core_->Drive();

            return true;
        }

        void set_parallel_logger(SP_ParallelLogger log) {
            this->log_ = log;
            this->driver_->set_paralog(log);
        }

        //--------------------------------------------------------------------------------------------------------------

        BraidCore* braid_core_ = nullptr;
        SP_ParallelLogger log_;
        SP_BraidBaseDriver driver_ = SPNULL;
        SP_SpaceTimeCommunicator comm_ = SPNULL;
        BraidSettings braid_settings_;

        //--------------------------------------------------------------------------------------------------------------
    };
}}


#endif