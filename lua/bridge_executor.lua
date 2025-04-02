
util = {}
util.xbraid = {}

function util.xbraid.create_braid_executor(desc,inst)
    method = BraidExecutor()
    if inst.initializer ~= nil then
        method:set_initializer(inst.initializer)
    end
    -- method:set_norm_provider() delegates to driver
    method:set_residual(desc.use_residual)
    -- level
    -- todo level    method:set_n_relax(level, desc.n_relax[level])
    -- todo level    method:set_c_factor(level, desc.c_factor[level])
    -- todo level    method:set_c_relax_weight(level, desc.c_relax_weight[level])
    method:set_max_levels(desc.max_levels)
    method:set_skip_downcycle_work(desc.skip_downcycle_work)
    method:set_min_coarse(desc.min_coarse)

    method:set_max_iterations(desc.conv.max_iterations)
    method:set_absolute_tol(desc.conv.absolut_tol)
    method:set_relative_tol(desc.conv.relative_tol)

    method:set_temporal_norm(desc.temporal_norm)
    method:set_sequential(desc.sequential)
    method:set_store_values(desc.store_values)

    method:set_spatial_coarsen_and_refine(desc.spatial_coarsen_and_refine)

    method:set_refine(desc.refine)
    method:set_max_refinements(desc.max_refinements)

    method:set_access_level(desc.access_level)
    method:set_print_level(desc.print_level)
    method:set_print_file(desc.print_file)
    if type(desc.cycle_type) == "string" then
        method:set_cycle_type(desc.cycle_type)
    end
    if desc.cycle.f_cycle then
        method:set_cycle_fmg()
        method:set_cycle_nfmg(desc.cycle.nfmg)
        method:set_cycle_nfmgv(desc.cycle.nfmgv)
    end

    if desc.sync then
        method:set_sync()
    end
    if desc.increase_max_levels then
        method:set_increase_max_levels()
    end

    method:set_relax_only_cg(desc.relax_only_cg)
    method:set_agg_c_factor(desc.agg_c_factor)
    method:set_periodic(desc.periodic)
    method:set_final_fc_relax(desc.final_fc_relax)
    method:set_reverted_ranks(desc.reverted_ranks)
    method:set_richardson_estimation(desc.richardson_estimation, desc.use_extrapolation, desc.local_order)
    method:set_file_io_level(desc.file_io_level)

    method:set_t_points_cutoff(desc.t_points_cutoff)
    if desc.full_residual_norm then
        method:set_full_residual_norm()
    end

    method:set_time_grid() -- todo empty method

    method:set_app(inst.driver) -- todo set app
    method:set_parallel_logger(inst.logger) -- todo set app
end