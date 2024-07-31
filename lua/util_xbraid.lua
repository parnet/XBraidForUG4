



-- --------------------------------------------------------------------------------
-- poro
-- --------------------------------------------------------------------------------

util = util or {}
util.xbraid = util.xbraid or {}
util.factory = util.factory or {}

print("Register XBraid-Lua Scripts")

function util.xbraid.cfactor_from_string(str_cfactor)
    print("lua  - cfactor_from_string")
    -- zero index is finest level
    cfactor = {}
    i = 1
    for s in str_cfactor:gmatch("[^_]+") do
        cfactor[i] = tonumber(s)
        i = i + 1
    end
    print()
    return cfactor, i
end

function util.xbraid.create_ThetaStepper(level, level_desc, inst)
    print("lua  - createThetaStepper")
    integrator = ThetaStepper()
    integrator:set_theta(level_desc.theta)
    integrator:set_reassemble_threshold(level_desc.threshold)
    integrator:set_domain(inst.domain)
    integrator:set_solver(lsolver)
    print()
    return integrator
end

function util.xbraid.set_approx_space(desc,inst)
    print("lua  - set_approx_space")
    print()
end

function util.xbraid.set_initializer(desc,inst,initializer)
    print("lua  - set_initializer")
    print()
end

function util.xbraid.set_relax_and_cycle(desc,inst)
    print("lua  - set_relax_and_cycle")
    if inst.braid == nil then
        print("    Create Braid Instance first!")
        exit()
    end
    short = ""
    cycle_type = desc.mgrit_cycle_type
    if cycle_type == "F" then
        print("    MGRIT uses F-Cycle")
        short = short .."F "
        inst.braid:set_cycle_fmg()
    elseif cycle_type == "V" then
        print("    MGRIT uses V-Cycle")
        short = short .. "V "
    else
        print("    Invalid MGRIT-cycle parameter - using V cycle")
        exit()
    end

    relax_type = desc.mgrit_relax_type
    if relax_type == "F" then
        inst.braid:set_n_relax(-1, 0)
        inst.braid:set_n_relax(0, 0)
        print("    MGRIT uses F-Relaxation on all level")
        short = short .. "F"
    elseif relax_type == "FCF" then
        inst.braid:set_n_relax(-1, 1)
        inst.braid:set_n_relax(0, 1)
        print("    MGRIT uses FCF-Relaxation on all level")
        short = short .. "FCF"
    elseif relax_type == "FFCF" then
        inst.braid:set_n_relax(-1, 1)
        inst.braid:set_n_relax(0, 0)
        short = short .. "F-FCF"
        print("    MGRIT uses F-FCF-Relaxation")
    elseif relax_type == "FCFF" then
        inst.braid:set_n_relax(-1, 0)
        inst.braid:set_n_relax(0, 1)
        short = short .. "FCF-F"
        print("    MGRIT uses FCF-Relaxation on coarse level and F-Relaxation on fines level ")
    elseif relax_type == "user" then
        if desc.nrelax_default ~= nil then
            inst.braid:set_n_relax(-1, desc.nrelax_default)
        else
            inst.braid:set_n_relax(-1, 0)
        end
        if desc.nrelax == "table" then
            for key,value in desc.nrelax do
                inst.braid:set_n_relax(key, value)
            end
            short = short .. " user-defined"
        end
    else
        print("Invalid MGRIT-relax parameter - using FCF relaxation")
        exit()
    end
    print()

end

function util.xbraid.write_script(desc,inst)
    print("lua  - write_script")
    if desc.write_script then
        paralog_script = Paralog()
        paralog_script:set_comm(inst.space_time_communicator)
        paralog_script:set_file_name("script")
        paralog_script:init()

        inst.braid:set_paralog_script(paralog_script)
    end
    print()
end

function util.xbraid.set_conv_check(desc, inst)
    print("lua  - set_conv_check")
    inst.braid:set_max_iterations(desc.conv_check.max_iter)
    if (desc.conv_check.absolute ~= 0) then
        inst.braid:set_absolute_tol(desc.conv_check.absolute)
    end

    if (desc.conv_check.reduction ~= 0) then
        inst.braid:set_relative_tol(desc.conv_check.reduction)
    end

    inst.prepared_conv_check = true
    print()
end

function util.xbraid.create_time_hierarchy(desc,inst)
    print("lua  - create_time_hierarchy")
    print("    Create time hierarchy!: ")
    if desc.default_cfactor ~= nil then
        inst.braid:set_c_factor(-1, desc.default_cfactor)
    end

    if type(desc.cfactor) == "table" then
        for key, value in pairs(desc.cfactor) do
            inst.braid:set_c_factor(key - 1, value)
        end
    else
        inst.braid:set_c_factor(-1, desc.cfactor)
    end

    max_level = desc.max_level

    fnum_time = desc.time.numtime
    print("    \tlvl  | \tcfac | \tnum-time")
    print("    \t--------------------------")
    base_reached = false
    for i = 1, #cfactor do
        cfac = cfactor[i]
        fnum_time_p = fnum_time / cfac
        if fnum_time_p < 2 or i == max_level then
            base_reached = true
            print("    \t",i-1, ": \tbase","\t", fnum_time)
            num_level = i
            break
        else
            print("    \t",i-1, ": \t", cfac,"\t",fnum_time)
        end
        fnum_time = fnum_time_p
        num_level = i
    end

    if not base_reached then
        print("    \t",i-1, ": \t", fnum_time, "\t base")
    end
    print()
    print("    Time hierarchy created: num level = ",num_level)
    inst.prepared_time_hierarchy = true
    inst.num_level = num_level
    print()

end

function util.factory.create_timedisc(desc)
    print("lua  - create_timedisc")
    if type(desc.integrator) == "string" then -- only default values
        if type(desc.conv_check) == "string" then -- only default values
        
        end
        if desc == "ThetaTimeStep" then
            --integrator:set_num_steps(num_steps) -> FixedStepThetaIntegrator
            --desc.solver == NL => NLFixedStepThetaIntegrator, NLThetaIntegrator
            -- orderOrTheta => BDF_Integrator
            integrator = ThetaSingleTimeStep()
            integrator:set_theta(1)
            integrator:set_reassemble_threshold()
            integrator:set_domain(inst.domain)
            integrator:set_solver(inst.lsolver)
        end
    end
    print()
end

--integrator:set_num_steps(num_steps) -> FixedStepThetaIntegrator
--desc.solver == NL => NLFixedStepThetaIntegrator, NLThetaIntegrator
-- orderOrTheta => BDF_Integrator

function util.xbraid.get_integrator_from_desc(desc,inst, level_desc)
    print("lua  - get_integrator_from_desc")
    if type(level_desc.integrator) == "userdata" then
        inst:set_integrator(desc.default.integrator) -- solver, conv_check and integrator included
    else
        solver = nil
        if integrator_desc.solver == "userdata" then
            solver = integrator_desc.solver
        else
            --todo get from description
        end

        if integrator_desc.solver == "userdata" then
            solver = integrator_desc.solver
        else
            --todo get from description
        end

        if integrator_desc.integrator == "ThetaTimeStep" then

        end
    end
    print()
end

function util.xbraid.get_solver(desc,inst,level_desc,level_inst)
    print("lua  - get_solver")
    if type(level_desc.solver) == "userdata" then
        return level_desc.solver
    else
        print("Not implementet yet!") --todo implement
    end
    exit()
end

function util.xbraid.get_conv_check(desc,inst,level_desc,level_inst)
    print("lua  - get_conv_check")
    if type(level_desc.conv_check) == "userdata" then
        return level_desc.conv_check
    else
        print("Not implementet yet!") --todo implement
        exit()
    end
end



function util.xbraid.get_integrator(desc,inst,level_desc, level_inst)
    print("lua  - get_integrator")
    if type(level_desc.integrator) == "userdata" then
        return level_desc.conv_check
    elseif  type(level_desc.integrator) == "string" then
        integrator = nil
        if level_desc.integrator == "ThetaTimeStep" then
            integrator = ThetaSingleTimeStep()
            integrator:set_domain(inst.domain_disc)
            integrator:set_theta(1)
            integrator:set_reassemble_threshold(1e-14)
            integrator:set_solver(level_inst.solver)

        elseif level_desc.integrator == "ThetaIntegrator" then
            integrator = ThetaIntegrator()
            integrator:set_domain(inst.domain_disc)
            integrator:set_theta(1)
            integrator:set_reassemble_threshold(1e-14)
            integrator:set_solver(level_inst.solver)

        elseif level_desc.integrator == "ThetaIntegratorNL" then
            integrator = ThetaIntegratorNL()
            integrator:set_domain(inst.domain_disc)
            integrator:set_order(1)
            integrator:set_solver(level_inst.solver)

        else
            print(level_desc.integrator)
            print("Integrator type must be specified by instance (userdata) or by a table")
            exit()
        end
        return integrator
    elseif type(level_desc.integrator) == "table" then
        print("Not implementet yet! (table) ") --todo implement
        exit()
    else
        print("Not implementet yet! (other) ") --todo implement
        exit()
    end
end

function util.xbraid.create_hierarchy_simple(desc,inst)
    print("lua  - create_hierarchy_simple")
    inst.driver = driver

    level_inst = {}
    level_inst.conv_check = util.xbraid.get_conv_check(desc,inst,desc.default,level_inst)
    level_inst.solver = util.xbraid.get_solver(desc,inst,desc.default,level_inst)
    level_inst.integrator = util.xbraid.get_integrator(desc,inst,desc.default,level_inst)

    print(desc.default.integrator)

    if desc.driver == "BasicDriver" then
        print("BasicDriver")
        driver:set_default_integrator(level_inst.integrator)

    elseif desc.driver == "NonlinearIntegrator" then
        print("BasicDriver")

    elseif desc.driver == "Integrator" then
        print("BasicDriver")
        driver:set_default_integrator()

    elseif desc.driver == "IntegratorFactory" then
        print("BasicDriver")
        driver:set_fine_time_integrator()
        driver:set_coarse_time_integrator()

    elseif desc.driver == "ResidualStepper" then
        print("ResidualStepper")
    else
        print("    Driver is unknown")
        exit()
    end
    -- get solver
    -- get conc_check
    -- get integrator
    -- create app
    -- set values
    print()
end

function util.xbraid.create_hierarchy(desc,inst)
    print("lua  - create_hierarchy")

    if desc.level_config == "simple" then
        print("    simple level configuration")
        util.xbraid.create_hierarchy_simple(desc,inst)
    end


    num_level = inst.num_level
    if type(desc.solver) == "table" then
        for i = 1,num_level do
            print(i)
            --if desc then
            --    app:set_integrator(k, xbraid_util.createThetaStepper(domainDiscT, lsolver, theta, threshold))
            --end
        end
    else
        -- app:set_default_integrator()
    end
    print("    Created Solver")
    print()
end

function util.xbraid.create_observer(observer, desc)
    print("lua  - create_observer")
    for i = 1, #desc.cfactor do
        observer:set_c_factor(i - 1, desc.cfactor[i])
    end
    print()
end

function util.xbraid.set_domain_disc(desc,inst,domain_disc)
    inst.domain_disc = domain_disc
    print()
end

function util.xbraid.create_driver(desc,inst)
    print("lua  - create_driver")
    if inst.domain_disc == nil then
        print("    use set_domain() to set the domain before calling this function")
        exit()
    end

    driver = nil
    if desc.driver == "BasicDriver" then
        driver = BasicDriver()
    elseif desc.driver == "NonlinearIntegrator" then
        driver = BraidNLIntegrator()
    elseif desc.driver == "Integrator" then
        driver = BraidIntegrator()
    elseif desc.driver == "IntegratorFactory" then
        driver = BraidIntegratorFactory()
    elseif desc.driver == "ResidualStepper" then
        driver = BraidResidualStepper()
    else
        print("    Driver is unknown")
        exit()
    end


    -- driver:set_verbose(desc.verbose) -- todo deleted
    driver:set_start_time(desc.time.t_0)
    driver:set_end_time(desc.time.t_end)
    driver:set_number_of_timesteps(desc.time.numtime)
    driver:set_time_values(desc.time.t_0, desc.time.t_end, desc.time.numtime)
    driver:set_max_levels(desc.max_level)
    -- driver:set_domain(domain)
    if desc.observer ~= nil then
        driver:set_observer(desc.observer)
    end

    if desc.xb_observer ~= nil then
        driver:set_xb_observer(desc.xb_observer)
    end

    inst.driver = driver
    print("    ",desc.driver, " created as driver for XBraid")
    print()
    return driver
end

function util.xbraid.adaptive(desc,inst) -- todo check and use ?
    print("lua  - adaptive")
    inst.driver:set_ref_factor(desc.rich_refine)
    inst.driver:set_threshold(desc.rich_bound)
    inst.driver:set_tol(1e-3, 1e-14)(desc.rich_bound)
    print()
end

function util.xbraid.set_communicator(desc,inst, spc)
    print("lua  - set_communicator")
    inst.communicator = spc
    print()
end

function util.xbraid.set_norm(desc,inst)
    print("lua  - set_norm")
    if desc.norm then
        braid:set_norm_provider(desc.norm)
    else
        print("    Norm provider must be used")
        exit()
    end
    print()
end

function util.xbraid.create_instance(desc,inst)
    print("lua  - create_instance")
    braid = BraidExecutor(inst.communicator, inst.driver)
    inst.braid = braid
    util.xbraid.set_relax_and_cycle(desc, inst) -- todo  <----
    braid:set_residual(desc.use_residual)
    print("    Residual: ", desc.use_residual)
    braid:set_temporal_norm(desc.temporal_norm)
    braid:set_access_level(desc.access_level)
    braid:set_print_level(desc.print_level)
    braid:set_store_values(desc.store_values)
    braid:set_skip_downcycle_work(desc.skip_downcycle_work)
    braid:set_max_levels(desc.max_level)
    braid:set_min_coarse(desc.min_coarsening)
    braid:set_sequential(desc.sequential)
    braid:set_spatial_coarsen_and_refine(desc.spatial_coarsen_and_refine)
    braid:set_refine(desc.time_refinement)
    braid:set_max_refinements(desc.max_refinement)
    braid:set_print_file(desc.printfile)
    -- braid:set_filename(desc.outputfile) -- todo deleted

    braid:set_richardson_estimation(desc.richardson_estimation,
                                    desc.richardson_extrapolation,
                                    desc.richardson_local_order)

    if inst.prepared_time_hierarchy == nil then
        util.xbraid.create_time_hierarchy(desc,inst)
    end
    if inst.prepared_conv_check == nil then
        util.xbraid.set_conv_check(desc, inst)
    end
    if desc.sync then
        braid:set_sync()
    end
    --inst.driver:init()
    --braid:print_settings()
    --inst.driver:print_settings()
    --inst.braid = braid
    util.xbraid.set_norm(desc,inst)
    util.xbraid.create_hierarchy(desc,inst)
    print()
    return braid
end


