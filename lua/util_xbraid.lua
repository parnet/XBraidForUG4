-- --------------------------------------------------------------------------------
-- poro
-- --------------------------------------------------------------------------------

util = util or {}
util.xbraid = util.xbraid or {}
util.factory = util.factory or {}

function util.xbraid.split_communicator()
    num_spatial_procs = util.GetParamNumber("--npx", 1, "number of spatial processors (must divide number of total_processors number)")
    threads = util.GetParamNumber("--threads", 1, "number of spatial processors (must divide number of total_processors number)")

    num_world_ranks = NumProcs()
    space_time_communicator = SpaceTimeCommunicator()

    if num_world_ranks % num_spatial_procs == 0 then
        space_time_communicator:split(num_spatial_procs)
        local num_temporal_procs = num_world_ranks / num_spatial_procs;
        local num_spatial_procs = num_spatial_procs
        print("temporal x spatial = world")
        print(num_temporal_procs.." x " .. num_spatial_procs .. " = " .. num_world_ranks)
    else
        space_time_communicator:split(1)
        print("temporal x spatial = world")
        print(num_world_ranks.." x " .. 1 .. " = " .. num_world_ranks)
    end
    -- space_time_communicator:set_openmp(0,threads)
    return space_time_communicator
end

function util.xbraid.cfactor_from_string(str_cfactor)
    local cfactor = {}
    local i = 1
    for s in str_cfactor:gmatch("[^_]+") do
        cfactor[i] = tonumber(s)
        i = i + 1
    end
    return cfactor, i
end

function util.xbraid.numref_from_string(str_numref)
    local numrefs = {}
    local i = 1
    for s in str_numref:gmatch("[^_]+") do
        numrefs[i] = tonumber(s)
        i = i + 1
    end
    return numrefs, i
end

function util.xbraid.create_ThetaStepper(level, level_desc, inst)
    integrator = ThetaStepper()
    integrator:set_theta(level_desc.theta)
    integrator:set_reassemble_threshold(level_desc.threshold)
    integrator:set_domain(inst.domain)
    integrator:set_solver(lsolver) -- (ø)[[todo missing argument?]]
    return integrator
end

function util.xbraid.set_approx_space(desc,inst)
    print("braid  - set_approx_space")
    print("(ø)[[todo - empty method]]")
end

function util.xbraid.set_initializer(desc,inst,initializer)
    if type(initializer) == "userdata" then
        inst.driver:set_initializer(initializer)
    elseif initializer == nil then
        print("Default Initializer will be used and set at a later point")
        return -- use default ( start value generator )
    else
        print("implementation missing")
        exit()
    end
end

function util.xbraid.set_relax_and_cycle(desc,inst)
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
    if desc.write_script then
        paralog_script = Paralog()
        paralog_script:set_comm(inst.space_time_communicator)
        paralog_script:set_file_name("script")
        paralog_script:init()

        inst.braid:set_paralog_script(paralog_script)
    end
end

function util.xbraid.set_conv_check(desc, inst)
    inst.braid:set_max_iterations(desc.conv_check.max_iter)
    if (desc.conv_check.absolute ~= 0) then
        inst.braid:set_absolute_tol(desc.conv_check.absolute)
    end
    if (desc.conv_check.reduction ~= 0) then
        inst.braid:set_relative_tol(desc.conv_check.reduction)
    end
    inst.prepared_conv_check = true
end

function util.xbraid.create_time_hierarchy(desc,inst)
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

    if type(desc.level_num_ref) == "table" then
        for key, value in pairs(desc.level_num_ref) do
            -- todo add for NL driver -- inst.driver:set_level_num_ref(key - 1, value) -- todo COARSE
        end
    end

    max_level = desc.max_level

    fnum_time = desc.time.numtime
    print("    \tlvl  | \tcfac | \tnum-time | \tnum-ref")
    print("    \t--------------------------------------")
    base_reached = false
    for i = 1, #desc.cfactor do
        cfac = desc.cfactor[i]
        numref = desc.level_num_ref[i]
        fnum_time_p = fnum_time / cfac
        if fnum_time_p < 2 or i == max_level then
            base_reached = true
            print("    \t",i-1, ": \tbase","\t", fnum_time,"\t\t",numref)
            num_level = i
            break
        else
            print("    \t",i-1, ": \t", cfac,"\t",fnum_time,"\t\t",numref)
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
    -- exit()
    -- ------------------------------
    local transfer = StdTransfer()
    transfer:enable_p1_lagrange_optimization(true)
    -- todo COARSE inst.driver:set_transfer(transfer)
    -- todo COARSE inst.driver:set_restrict(RestrictWrapper())
    -- todo COARSE inst.driver:set_prolongate(ProlongateWrapper())
    -- ------------------------------
end


function util.xbraid.get_solver(desc,inst,level_desc,level_inst)
    print("braid  - get_solver")
    print("level_desc.solver=",type(level_desc.solver))
    print("solver type",type(level_desc.solver))
    if type(level_desc.solver) == "userdata" then
        return level_desc.solver
    else
        print(type(level_desc))
        print("----> util.xbraid.get_solver ::: Not implementet yet!")
        print()
        print()
        exit()
    end
end

--- todo describe what function does
--- @desc
--- @inst
--- @level_desc
--- @level_inst
function util.xbraid.get_conv_check(desc,inst,level_desc,level_inst)
    print("braid - util.xbraid.get_conv_check")
    if type(level_desc.conv_check) == "userdata" then
        return level_desc.conv_check
    else
        print("util.xbraid.get_conv_check ::: Not implementet yet!") --todo implement
        exit()
    end
end



function util.xbraid.get_integrator(desc,inst,level_desc, level_inst)
    print("braid  - get_integrator")
    if type(level_desc.integrator) == "userdata" then
        return level_desc.conv_check
    elseif  type(level_desc.integrator) == "string" then
        integrator = nil
        if level_desc.integrator == "ThetaTimeStep" then
            print("Time Integrator: ThetaTimeStep")
            integrator = ThetaSingleTimeStep()
            integrator:set_domain(inst.domain_disc)
            integrator:set_theta(1)
            integrator:set_reassemble_threshold(1e-14)
            integrator:set_solver(level_inst.solver)

        elseif level_desc.integrator == "ThetaIntegrator" then
            print("Time Integrator: ThetaIntegrator")
            integrator = ThetaIntegrator()
            integrator:set_domain(inst.domain_disc)
            integrator:set_theta(1)
            integrator:set_reassemble_threshold(1e-14)
            integrator:set_solver(level_inst.solver)

        elseif level_desc.integrator == "ThetaIntegratorNL" then
            print("Time Integrator: ThetaIntegrator Nonlinear")
            integrator = ThetaIntegratorNL()
            integrator:set_domain(inst.domain_disc)
            integrator:set_theta(1)
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
    print("braid  - create_hierarchy_simple")
    inst.driver = driver

    print(desc)

    level_inst = {}
    level_inst.conv_check = util.xbraid.get_conv_check(desc,inst,desc.default,level_inst)
    level_inst.solver = util.xbraid.get_solver(desc,inst,desc.default,level_inst)
    level_inst.integrator = util.xbraid.get_integrator(desc,inst,desc.default,level_inst)

    print(desc.default.integrator)

    if desc.driver == "BasicDriver" then
        print("Driver: Basic Driver")
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
    print()
end

function util.xbraid.create_hierarchy_finecoarse(desc,inst)
    print("braid  - create_hierarchy_fine_coarse")
    inst.driver = driver

    print(desc.fine)
    print(desc.fine.solver)
    print(desc.fine.conv_check)
    print(desc.fine.integrator)
    print()
    print()
    print()
    print()
    fine_level_inst = {}
    fine_level_inst.conv_check = util.xbraid.get_conv_check(desc,inst,desc.fine,fine_level_inst)
    fine_level_inst.solver = util.xbraid.get_solver(desc,inst,desc.fine,fine_level_inst)
    fine_level_inst.integrator = util.xbraid.get_integrator(desc,inst,desc.fine,fine_level_inst)
    driver:set_integrator(0,fine_level_inst.integrator)

    for level = 1, desc.max_level - 1 do
        coarse_level_inst = {}
        coarse_level_inst.conv_check = util.xbraid.get_conv_check(desc,inst,desc.coarse,coarse_level_inst)
        coarse_level_inst.solver = util.xbraid.get_solver(desc,inst,desc.coarse,coarse_level_inst)
        coarse_level_inst.integrator = util.xbraid.get_integrator(desc,inst,desc.coarse,coarse_level_inst)
        driver:set_integrator(level,coarse_level_inst.integrator)
        print("added to level lvl = " .. level)
    end

    print()
end


function util.xbraid.create_hierarchy(desc,inst)
    print("braid  - create_hierarchy")

    if desc.level_config == "simple" then
        util.xbraid.create_hierarchy_simple(desc,inst)
    elseif desc.level_config == "finecoarse" then
        util.xbraid.create_hierarchy_finecoarse(desc,inst)
    elseif desc.level_config == "leveldependend" then
        util.xbraid.create_hierarchy_leveldependend(desc,inst)
    else
        print("invalid argument for level configuration [ simple | finecoarse | leveldependend ]")
    end

end

function util.xbraid.create_observer(desc, inst, vtk, scriptor)
    -- print("braid  - create_observer")
    -- for i = 1, #desc.cfactor do
    --     observer:set_c_factor(i - 1, desc.cfactor[i])
    -- end
    -- print()
    --vtk_observer = VTK_Observer(vtk,"solution")
    --inst.driver:attach_observer(vtk_observer)

    --vtk_processobserver = VTK_ProcessObserver(vtk,"process")
    inst.driver:attach_xbraid_observer(scriptor)
    -- exit()
end

function util.xbraid.set_domain_disc(desc,inst,domain_disc)
    if domain_disc == nil then
        print("util.xbraid.set_domain_disc ::: Domain Disc is nil")
        exit()
    end
    inst.domain_disc = domain_disc
    print()
end

function util.xbraid.create_driver(desc,inst)
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
        print("    Driver argument is unknown")
        exit()
    end


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
    print("braid  - adaptive")
    inst.driver:set_ref_factor(desc.rich_refine)
    inst.driver:set_threshold(desc.rich_bound)
    inst.driver:set_tol(1e-3, 1e-14)(desc.rich_bound)
    print()
end

function util.xbraid.set_communicator(desc, inst, spc)
    inst.communicator = spc
end

function util.xbraid.set_norm(desc,inst)
    if type(desc.norm) == "userdata" then
        inst.driver:set_norm_provider(desc.norm)
    elseif type(desc.norm) == "string" then
        if desc.norm == "L2" then
            norm = BraidEuclidianNorm();
            inst.driver:set_norm_provider(norm)
        else
            print("util.xbraid.set_norm ::: ", desc.norm, " currently not implemented ")
        end
    else
        print("    Norm provider must be used")
        exit()
    end
    print()
end

function util.xbraid.create_instance(desc,inst)
    braid = BraidExecutor(inst.communicator, inst.driver)
    inst.braid = braid
    util.xbraid.set_relax_and_cycle(desc, inst) -- todo  <----
    braid:set_residual(desc.use_residual)
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
    util.xbraid.set_norm(desc, inst)
    util.xbraid.create_hierarchy(desc,inst)
    print()
    return braid
end


