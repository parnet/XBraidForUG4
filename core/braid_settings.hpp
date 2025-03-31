#ifndef UGPLUGIN_XBRAIDFORUG4_CORE_BRAID_SETTINGS_HPP
#define UGPLUGIN_XBRAIDFORUG4_CORE_BRAID_SETTINGS_HPP


namespace ug{ namespace xbraid {

        class BraidSettings {
        public:
            //----------------------------------------------------------------------------------------------------------

            bool m_increase_max_level = false;
            bool m_skip = false;
            bool m_cycle_fmg = false;
            bool m_residual = false;
            bool m_sync = false;
            bool m_final_fc_relax = false;
            bool finished = false;
            bool m_file_io_level = true;
            bool m_relax_only_cg = false;
            bool m_est_error = false; // richardson-based error estimation
            bool m_richardson = false; //  richardson-based extrapolation for finest grid
            bool m_sequential = false;
            bool m_full_r_norm = false;
            bool m_coarsen_and_refine = false;
            bool m_refine = false;
            int eval_level = 0;
            int m_max_level = 30;
            int m_min_coarse = -1;
            int m_max_iter = 100;
            int m_temp_norm = 2;
            int m_store_values = -1;
            int m_t_points_cutoff = -1;
            int m_max_refinements = 200;
            int m_access_level = 1;
            int m_print_level = 2;
            int m_cycle_nfmg = -1;
            int m_cycle_nfmgv = 1;
            int m_reverted_ranks = 0;
            int m_local_order = -1; // 2 for bachward euler
            int m_periodic = -1;
            double m_abs_tol = -1.0;
            double m_rel_tol = -1.0;
            const char *m_print_file = "braid_runtime.out";
            int m_n_relax_default = 1;
            std::map<int, int> m_n_relax;
            int m_c_factor_default = 2;
            std::map<int, int> m_c_factor;
            double m_c_relax_weight_default = 1;
            std::map<int, double> m_c_relax_weight;

            //----------------------------------------------------------------------------------------------------------

            BraidSettings() = default;
            ~BraidSettings() = default;

            //----------------------------------------------------------------------------------------------------------
        };
    }
}
#endif