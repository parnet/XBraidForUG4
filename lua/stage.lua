inst = {
    domain_disc  = nil,
    approx_space = nil,
    time_disc = nil,

    initializer = nil, -- optional
    driver = 0,
    logger = 0,

    linear_solver = 0,
    nonlinear_solver = 0,

    transfer = 0,
    restriction = 0,
    prolongation = 0,


    observer = 0,
    vtk = 0,

    error_estimator = 0,
    comm = nil,


    level = {[3] = {linear_solver = 0,
                    nonlinear_solver = 0}}
}
-- transfer or (restriction and prolongation) needed


desc = {
    time_interval = {
        start_time,
        end_time,
        time_steps
    },
    use_residual,
    hierarchy = {
        max_levels,
        store_values,
    },

    min_coarse,
    level_config = "simple", -- simple, finecoarse, leveldependent
    temporal_norm = 3 , -- 1 sum abs(x_i) ,2 sqrt(sum x_i**2), 3 max(abs(x_i))

    --cat:other
    sequential,
    --value::verbose
    spatial_coarsen_and_refine,
    max_refinements,
    temporal_refine = { -- or error estimate
        richardson_refine,
        richarson_threshold
    },
    richardson_estimation,
    use_extrapolation,

    refine,

    -- file settings
    access_level,
    print_level,
    print_file,
    file_io_level,
    --

    spatial_coarsen = {
        lagrange_optimization = true
    },

    sync,
    increase_max_levels,

    agg_c_factor,
    periodic,
    final_fc_relax,
    reverted_ranks,


    local_order,

    t_points_cutoff,
    full_residual_norm,

    -- cycle
    -- todo search and replace  cycle_type   ->    cycle.cycle_type , -- V  or F
    skip_downcycle_work,
    -- relaxation
    -- crelaxwt (default = 1.0)
    relax_only_cg,
    relaxation = {
        type = "F", -- "F", "FCF", "FFCF", "FCFF", "user"
        default = 0, -- 0=F Relaxation, 1=FCF Relaxation
        level = nil -- or table {1: 4, 2: 8 etc} with additional CF-Relaxations per level -- todo from string ? 1_2_2_2
    },
    cycle = {
        relax_type = "F", -- "F" or "FCF"
        cycle_type = "V", -- "V" or "F"
        nfmg = nil , -- number of f-cycles before switching to V cycles
        nfmgv = 1, -- number of v-cycles at each level (default = 1)
    },


    integrator = {
        type = "method", -- "method" or "factory"
        num_steps = 1,  -- for ConstStepLinearTimeIntegratorFactory
        dt_max = 1e-2,  -- for LimexFactory, SimpleIntegratorFactory
        dt_min = 1e-10, -- for LimexFactory, SimpleIntegratorFactory
        tol = 1e-3,     -- for LimexFactory
    },

    conv_check = {
        max_iterations,
        absolut_tol,
        relative_tol
    },
    newton_conv_check = { -- class --> ConvCheck & NewtonSolver
        maximum_steps = 1,
        reduction = 0.9999999999999999999,
        minimum_defect = 1e-14,
        verbose = true,
        suppress_unsuccessful = true
    },
    linear_conv_check = {
        maximum_steps = 100,
        reduction = 1e-8,
        minimum_defect = 1e-14,
        verbose = true,
        suppress_unsuccessful = true
    }
}