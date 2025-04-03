#ifndef UGPLUGIN_XBRAIDFORUG4_CORE_BRAID_SETTINGS_HPP
#define UGPLUGIN_XBRAIDFORUG4_CORE_BRAID_SETTINGS_HPP


namespace ug{ namespace xbraid {

        class BraidSettings {
        public:
            //----------------------------------------------------------------------------------------------------------

            bool increase_max_level_ = false;
            bool skip_ = false;
            bool cycle_fmg_ = false;
            bool residual_ = false;
            bool sync_ = false;
            bool final_fc_relax_ = false;
            bool finished_ = false;
            bool file_io_level_ = true;
            bool relax_only_cg_ = false;
            bool est_error_ = false; // richardson-based error estimation
            bool richardson_ = false; //  richardson-based extrapolation for finest grid
            bool sequential_ = false;
            bool full_r_norm_ = false;
            bool coarsen_and_refine_ = false;
            bool refine_ = false;
            int eval_level_ = 0;
            int max_level_ = 30;
            int min_coarse_ = -1;
            int max_iter_ = 100;
            int temp_norm_ = 2;
            int store_values_ = -1;
            int t_points_cutoff_ = -1;
            int max_refinements_ = 200;
            int access_level_ = 1;
            int print_level_ = 2;
            int cycle_nfmg_ = -1;
            int cycle_nfmgv_ = 1;
            int reverted_ranks_ = 0;
            int local_order_ = -1; // 2 for bachward euler
            int periodic_ = -1;
            int c_factor_default_ = 2;
            double abs_tol_ = -1.0;
            double rel_tol_ = -1.0;
            double c_relax_weight_default_ = 1;
            const char * print_file_ = "braid_runtime.out";
            int n_relax_default_ = 1;
            std::map<int, int> n_relax_;
            std::map<int, int> c_factor_;
            std::map<int, double> c_relax_weight_;

            //----------------------------------------------------------------------------------------------------------

            BraidSettings() = default;
            ~BraidSettings() = default;

            //----------------------------------------------------------------------------------------------------------
        };
    }
}
#endif