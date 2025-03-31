#ifndef UUGPLUGIN_XBRAIDFORUG4_CORE_BRAID_EXECUTOR_H
#define UGPLUGIN_XBRAIDFORUG4_CORE_BRAID_EXECUTOR_H

#include "bindings/lua/lua_user_data.h"
#include "driver/basic_driver.hpp"
#include "interface/observer_xbraid.hpp"
#include "space_time_communicator.hpp"
#include "braid_settings.hpp"
#include "initializer/start_value_initializer.hpp"

#include "util/pragma.hpp"



namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class BraidExecutor {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        using T_BraidBaseApp = BraidGridFunctionBase<TDomain, TAlgebra> ;
        using SP_BraidBaseApp = SmartPtr<T_BraidBaseApp> ;

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

        BraidCore* m_braid_core = nullptr; // todo change to SmartPtr?
        SP_ParallelLogger m_log;
        SP_BraidBaseApp m_driver = SPNULL;
        SP_SpaceTimeCommunicator m_comm = SPNULL;
        BraidSettings m_braid_settings;

        //--------------------------------------------------------------------------------------------------------------

        BraidExecutor(SP_SpaceTimeCommunicator& p_comm, SP_BraidBaseApp& p_app) {
            this->m_comm = p_comm;
            this->set_app(p_app);
        }

        virtual ~BraidExecutor() {
            delete m_braid_core;
        }

        //--------------------------------------------------------------------------------------------------------------

        void set_residual(bool residual) {
            if (residual ) { // && this->m_driver->can_residual_method
                this->m_braid_settings.m_residual = residual;
                this->m_braid_core->SetResidual();
            }
        }

        void set_n_relax(int level, int number) {
            if (level == -1) {
                this->m_braid_settings.m_n_relax_default = number;
            }
            else {
                this->m_braid_settings.m_n_relax[level] = number;
                if (level > this->m_braid_settings.eval_level) {
                    this->m_braid_settings.eval_level = level;
                }
            }
            this->m_braid_core->SetNRelax(level, number);
        }

        void set_c_factor(int level, int factor) {
            if (level == -1) {
                this->m_braid_settings.m_c_factor_default = factor;
            }
            else {
                this->m_braid_settings.m_c_factor[level] = factor;
                if (level > this->m_braid_settings.eval_level) {
                    this->m_braid_settings.eval_level = level;
                }
            }
            this->m_braid_core->SetCFactor(level, factor);
        }

        void set_max_levels(int maxLevel) {
            this->m_braid_settings.m_max_level = maxLevel;
            this->m_braid_core->SetMaxLevels(maxLevel);
            this->m_driver->set_max_levels(maxLevel);
        }

        void set_skip_downcycle_work(bool skip) {
            this->m_braid_settings.m_skip = skip;
            this->m_braid_core->SetSkip(skip);
        }

        void set_min_coarse(int minCoarse) {
            this->m_braid_settings.m_min_coarse = minCoarse;
            this->m_braid_core->SetMinCoarse(minCoarse);
        }

        void set_max_iterations(int max_iter) {
            this->m_braid_settings.m_max_iter = max_iter;
            this->m_braid_core->SetMaxIter(max_iter);
        }

        void set_absolute_tol(double tol) {
            this->m_braid_settings.m_abs_tol = tol;
            this->m_braid_core->SetAbsTol(tol);
        }

        void set_relative_tol(double tol) {
            this->m_braid_settings.m_rel_tol = tol;
            this->m_braid_core->SetRelTol(tol);
        }

        void set_temporal_norm(int nrm) {
            this->m_braid_settings.m_temp_norm = nrm;
            this->m_braid_core->SetTemporalNorm(nrm);
        }

        void set_sequential(bool sequential) {
            this->m_braid_settings.m_sequential = sequential;
            if (sequential) {
                this->m_braid_core->SetSeqSoln(1);
            } else {
                this->m_braid_core->SetSeqSoln(0);
            }
        }

        void set_store_values(int level) {
            this->m_braid_settings.m_store_values = level;
            this->m_braid_core->SetStorage(level);
        }

        void set_spatial_coarsen_and_refine(bool cnr) {
            if (cnr) {
                this->m_braid_settings.m_coarsen_and_refine = cnr;
                this->m_braid_core->SetSpatialCoarsenAndRefine();
            }
        }

        void set_refine(bool ref) {
            this->m_braid_settings.m_refine = ref;
            if (ref) {
                this->m_braid_core->SetRefine(1);
            } else {
                this->m_braid_core->SetRefine(0);
            }
        }

        void set_max_refinements(int number) {
            this->m_braid_settings.m_max_refinements = number;
            this->m_braid_core->SetMaxRefinements(number);
        }

        void set_access_level(int level) {
            this->m_braid_settings.m_access_level = level;
            this->m_braid_core->SetAccessLevel(level);
        }

        void set_print_level(int level) {
            this->m_braid_settings.m_print_level = level;
            this->m_braid_core->SetPrintLevel(level);
        }

        void set_default_print_file() {
            const char* file = "braid_runtime.out";
            this->set_print_file(file);
        }

        void set_print_file(const char* file) {
            this->m_braid_settings.m_print_file = file;
            this->m_braid_core->SetPrintFile(file);
        }

        void set_cycle_fmg() {
            this->m_braid_settings.m_cycle_fmg = true;
            this->m_braid_core->SetFMG();
        }

        void set_cycle_nfmg(int vu) {
            this->m_braid_settings.m_cycle_nfmg = vu;
            this->m_braid_core->SetNFMG(vu);
        }

        void set_cycle_nfmgv(int mu) {
            this->m_braid_settings.m_cycle_nfmgv = mu;
            this->m_braid_core->SetNFMGVcyc(mu);
        }

        void set_cycle_type(const char* ctyp) {
            std::cout << "debBraidExecutor::setCycleType(args)" << std::endl << std::flush;
            if (strcmp(ctyp, "V FCF") == 0) {
                this->m_braid_core->SetNRelax(-1, 1);
                this->m_braid_core->SetNRelax(0, 1);
            } else if (strcmp(ctyp, "V F") == 0) {
                this->m_braid_core->SetNRelax(-1, 0);
                this->m_braid_core->SetNRelax(0, 0);
            } else if (strcmp(ctyp, "V F-FCF") == 0) {
                this->m_braid_core->SetNRelax(-1, 1);
                this->m_braid_core->SetNRelax(0, 0);
            } else if (strcmp(ctyp, "F FCF") == 0) {
                this->m_braid_core->SetNRelax(-1, 1);
                this->m_braid_core->SetNRelax(0, 1);
                this->m_braid_settings.m_cycle_fmg = true;
                this->m_braid_core->SetFMG();
            } else if (strcmp(ctyp, "F F") == 0) {
                this->m_braid_core->SetNRelax(-1, 0);
                this->m_braid_core->SetNRelax(0, 0);
                this->m_braid_settings.m_cycle_fmg = true;
                this->m_braid_core->SetFMG();
            } else if (strcmp(ctyp, "F F-FCF") == 0) {
                this->m_braid_core->SetNRelax(-1, 1);
                this->m_braid_core->SetNRelax(0, 0);
                this->m_braid_settings.m_cycle_fmg = true;
                this->m_braid_core->SetFMG();
            } else {
                std::cout << "invalid cycle type" << std::endl;
            }
        }

        void set_sync() {
            this->m_braid_settings.m_sync = true;
            this->m_braid_core->SetSync();
        }

        void set_increase_max_levels() {
            this->m_braid_settings.m_increase_max_level = true;
            this->m_braid_core->SetIncrMaxLevels();
        }

        void set_relax_only_cg(bool setting) {
            this->m_braid_settings.m_relax_only_cg = setting;
            if (setting) {
                this->m_braid_core->SetRelaxOnlyCG(1);
            } else {
                this->m_braid_core->SetRelaxOnlyCG(0);
            }
        }

        void set_agg_c_factor(int cfactor0) {
            this->m_braid_core->SetAggCFactor(cfactor0);
        }

        void set_periodic(int periodic) {
            this->m_braid_settings.m_periodic = periodic;
            this->m_braid_core->SetPeriodic(periodic);
        }

        void set_final_fc_relax() {
            this->m_braid_settings.m_final_fc_relax = true;
            this->m_braid_core->SetFinalFCRelax();
        }

        void set_reverted_ranks(int ranks) {
            this->m_braid_settings.m_reverted_ranks = ranks;
            this->m_braid_core->SetRevertedRanks(ranks);
        }

        void set_richardson_estimation(bool use_richardson, bool use_extrapolation, int local_order) {
            this->m_braid_settings.m_est_error = use_richardson;
            this->m_braid_settings.m_richardson = use_extrapolation;
            this->m_braid_settings.m_local_order = local_order;
            this->m_braid_core->SetRichardsonEstimation(use_richardson, use_extrapolation, local_order);
        }

        void set_file_io_level(int level) {
            this->m_braid_settings.m_file_io_level = level;
            this->m_braid_core->SetFileIOLevel(level);
        }

        void set_c_relax_weight(int level, double weight) {
            if (level == -1) {
                this->m_braid_settings.m_c_relax_weight_default = weight;
            } else {
                this->m_braid_settings.m_c_relax_weight[level] = weight;
                if (level > this->m_braid_settings.eval_level) {
                    this->m_braid_settings.eval_level = level;
                }
            }
            this->m_braid_core->SetCRelaxWt(level, weight);
        }

        void set_t_points_cutoff(int cutoff) {
            this->m_braid_settings.m_t_points_cutoff = cutoff;
            this->m_braid_core->SetTPointsCutoff(cutoff);
        }

        void set_full_residual_norm() {
            this->m_braid_settings.m_full_r_norm = true;
            this->m_braid_core->SetFullRNormRes(_BraidAppResidual);
        }

        void set_time_grid() {} // todo implement or delete

        int get_num_iteration() {
            int iter = 0;
            this->m_braid_core->GetNumIter(&iter);
            return iter;
        }

        void get_c_factor() {} // todo implement or delete

        void get_residual_norms() {}  // todo implement or delete

        int get_num_level() {
            int number_of_level = 0;
            this->m_braid_core->GetNLevels(&number_of_level);
            return number_of_level;
        }

        int get_warm_restart() const {
            return this->m_braid_core->GetWarmRestart();
        }

        int get_distribution_lower() const {
            int lower;
            int upper;
            this->m_braid_core->GetDistribution(&lower, &upper);
            return lower;
        }

        int get_distribution_upper() const {
            int lower;
            int upper;
            this->m_braid_core->GetDistribution(&lower, &upper);
            return upper;
        }

        void get_distribution(int& lower, int& upper) const {
            this->m_braid_core->GetDistribution(&lower, &upper);
        }

        int get_id() const {
            int id = 0;
            this->m_braid_core->GetMyID(&id);
            return id;
        }

        braid_Core get_real_braid_core() const {
            return this->m_braid_core->GetCore();
        }

        BraidCore* get_braid_core() const {
            return this->m_braid_core;
        }

        void reset_core() {
            delete m_braid_core;
            m_braid_core = new BraidCore(this->m_comm->GLOBAL, this->m_driver.get()); // todo make smartptr?
        }

        void set_app(SP_BraidBaseApp p_app) {
            this->m_driver = p_app;
            this->m_driver->m_comm = this->m_comm;
            this->m_driver->comm_t = this->m_comm->TEMPORAL;
            this->reset_core();
        }

        SP_BraidBaseApp get_driver() {
            return this->m_driver;
        }

        void print_settings() {
            this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
            this->m_log->o << "Module Name:                                             Braid Executor" << std::endl;
            this->m_log->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----" << std::endl;
            this->m_log->o << "Max Level: " << this->m_braid_settings.m_max_level << std::endl;
            this->m_log->o << "Increase max level: " << this->m_braid_settings.m_increase_max_level << std::endl;
            this->m_log->o << "Skip Down Cylce Work: " << this->m_braid_settings.m_skip << std::endl;
            this->m_log->o << "Absolute Tol: " << this->m_braid_settings.m_abs_tol << std::endl;
            this->m_log->o << "Relative Tol: " << this->m_braid_settings.m_rel_tol << std::endl;
            this->m_log->o << "Temporal Norm: " << this->m_braid_settings.m_temp_norm << std::endl;
            this->m_log->o << "Max Iteration: " << this->m_braid_settings.m_max_iter << std::endl;
            this->m_log->o << "Sync: " << this->m_braid_settings.m_sync << std::endl;
            this->m_log->o << "Residual: " << this->m_braid_settings.m_residual << std::endl;
            this->m_log->o << "Cycle FMG: " << this->m_braid_settings.m_cycle_fmg << std::endl;
            this->m_log->o << "Cycle NFMG: " << this->m_braid_settings.m_cycle_nfmg << std::endl;
            this->m_log->o << "Cycle NFMGV: " << this->m_braid_settings.m_cycle_nfmgv << std::endl;
            this->m_log->o << "Store Values: " << this->m_braid_settings.m_store_values << std::endl;
            this->m_log->o << "Print File: " << this->m_braid_settings.m_print_file << std::endl;
            this->m_log->o << "Print Level: " << this->m_braid_settings.m_print_level << std::endl;
            this->m_log->o << "File IO Level: " << this->m_braid_settings.m_file_io_level << std::endl;
            this->m_log->o << "Access Level: " << this->m_braid_settings.m_access_level << std::endl;
            this->m_log->o << "Sequential: " << this->m_braid_settings.m_sequential << std::endl;
            this->m_log->o << "Coarsen & Refine: " << this->m_braid_settings.m_coarsen_and_refine << std::endl;
            this->m_log->o << "Min Coarse: " << this->m_braid_settings.m_min_coarse << std::endl;
            this->m_log->o << "Final FC Relax: " << this->m_braid_settings.m_final_fc_relax << std::endl;
            this->m_log->o << "Refine: " << this->m_braid_settings.m_refine << std::endl;
            this->m_log->o << "Max Refinements: " << this->m_braid_settings.m_max_refinements << std::endl;
            this->m_log->o << "Relax only CG: " << this->m_braid_settings.m_relax_only_cg << std::endl;
            this->m_log->o << "t-point cutoff: " << this->m_braid_settings.m_t_points_cutoff << std::endl;
            this->m_log->o << "Braid ID: " << this->get_id() << std::endl;
            this->m_log->o << "Distribution: " << this->get_distribution_lower() << " , " << this->get_distribution_upper() << std::endl;
            this->m_log->o << "Number of Level: " << this->get_num_level() << std::endl;
            this->m_log->o << "Number of Iteration: " << this->get_num_iteration() << std::endl;
            this->m_log->o << "Richardson-based Error Estimation: " << this->m_braid_settings.m_est_error << std::endl;
            this->m_log->o << "Richardson-based Extrapolation: " << this->m_braid_settings.m_richardson << std::endl;
            this->m_log->o << "Local order: " << this->m_braid_settings.m_local_order << std::endl;
            this->m_log->o << "Reverted Ranks: " << this->m_braid_settings.m_reverted_ranks << std::endl;
            this->m_log->o << "Periodic: " << this->m_braid_settings.m_periodic << std::endl;
            this->m_log->o << "Full residual norm: " << this->m_braid_settings.m_full_r_norm << std::endl;
            this->m_log->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----" << std::endl;
            this->m_log->o << "default" << m_braid_settings.m_n_relax_default << "    " << m_braid_settings.m_c_relax_weight_default << "    " << m_braid_settings.m_c_factor_default << "    " << std::endl;
            for (int i = 0; i <= this->m_braid_settings.eval_level; i++) {
                this->m_log->o << i << "    ";
                if (this->m_braid_settings.m_n_relax.find(i) != this->m_braid_settings.m_n_relax.end()) {
                    this->m_log->o << " " << this->m_braid_settings.m_n_relax[i] << " ";
                } else {
                    this->m_log->o << "(" << m_braid_settings.m_n_relax_default << ")";
                }
                this->m_log->o << "    ";
                if (this->m_braid_settings.m_c_relax_weight.find(i) !=
                    this->m_braid_settings.m_c_relax_weight.end()) {
                    this->m_log->o << " " << this->m_braid_settings.m_c_relax_weight[i] << " ";
                } else {
                    this->m_log->o << "(" << m_braid_settings.m_c_relax_weight_default << ")";
                }
                this->m_log->o << "    ";
                if (this->m_braid_settings.m_c_factor.find(i) != this->m_braid_settings.m_c_factor.end()) {
                    this->m_log->o << " " << this->m_braid_settings.m_c_factor[i] << " ";
                } else {
                    this->m_log->o << "(" << m_braid_settings.m_c_factor_default << ")";
                }
                this->m_log->o << std::endl;
            }
            this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        }

        void print_summary() {
            this->m_log->o << "----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----"
                << std::endl;
            this->m_log->o << "Number of Level: " << this->get_num_level() << std::endl;
            this->m_log->o << "Number of Iteration: " << this->get_num_iteration() << std::endl;
            int request = this->m_braid_settings.m_max_iter;
            auto* norms = new double[this->m_braid_settings.m_max_iter];
            this->m_braid_core->GetRNorms(&request, norms);
            for (int i = 0; i < request; i++) {
                this->m_log->o << i << ": " << norms[i] << std::endl;
            }
            this->m_log->o << "===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====" << std::endl;
        }

        void set_observer(SP_T_ITimeIntegratorObserver output) {
            this->m_driver->set_observer(output);
        }

        void set_xbraid_observer(SP_XBTIObserver output) {
            this->m_driver->set_xbraid_observer(output);
        }

        void set_initializer(SP_BraidInitializer initializer) {
            this->m_driver->m_initializer = initializer;
        }

        void set_norm_provider(SPSpatialNorm norm) {
            this->m_driver->m_norm = norm;
        }

        bool apply(SP_GridFunction u_stop, number t_stop, SP_GridFunction u0, number t0) {
            __debug(std::cout << "BraidExecutor::apply" << std::endl);

            this->m_driver->set_start_time(t0);
            this->m_driver->set_end_time(t_stop);
            this->m_driver->set_start_vector(u0);

            if (this->m_driver->m_initializer == SPNULL) {
                std::cout << "Create default Initializer " << std::endl;
                this->m_driver->m_initializer = make_sp(new GridFunctionInitializer<TDomain,TAlgebra>());
            }

            if (this->m_log == SPNULL) {
                auto logger = make_sp(new T_ParallelLogger());
                logger->set_comm(this->m_comm);
                logger->set_filename("std_output");
                logger->init();
                this->set_parallel_logger(logger);
            }

            this->m_driver->m_initializer->set_start_values(u0, t0);

            this->m_driver->restimate = this->m_braid_settings.m_est_error; // todo why here?
            this->m_driver->refine_time = this->m_braid_settings.m_refine; // todo why here?
            this->m_driver->init();
            this->m_braid_core->Drive();

            return true;
        }

        void set_parallel_logger(SP_ParallelLogger log) {
            this->m_log = log;
            this->m_driver->set_paralog(log);
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}


#endif