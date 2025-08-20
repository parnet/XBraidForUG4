-- EnvironmentParameter
np_x
np_t
np_total = npx * npt
redirect_output
environment
(method) method = "",

refine = {
        trefine
        time_refine
        time_max_refine
        t_points_cutoff
        }
--domain_disc
num_pre_refines
num_refs
-- spatial
spatial_refines

-- time hierarchy
set_time_grid
cfactor
        stringparameter 2_2_2_2_2
        listparameter {2,2,2,2,2}
cfactor_defaul
store_level store values on every coarse level up to value

-- convergence

spatial_norm l2, max, ...


-- basics

setbufallocfree
theta
order
gridstep
coarse max number of iteration for coarse solver
solver id
solver id coarse
level config


-- error estimation
rich_bound
rich_local_order
-- integrator
theta,
order,
grid_step,

-- other
coarse
min_coarsen_and_refine
min_coasening
level_config
level_numref
-- grid_name
-- grid_num_refs