




















function util.xbraid.set_relax_and_cycle(desc,inst)

    relax_type = desc.mgrit_relax_type
    if relax_type == "F" then
        inst.braid:set_n_relax(-1, 0)
        inst.braid:set_n_relax(0, 0)
        print("    MGRIT uses F-Relaxation on all level")
    elseif relax_type == "FCF" then
        inst.braid:set_n_relax(-1, 1)
        inst.braid:set_n_relax(0, 1)
        print("    MGRIT uses FCF-Relaxation on all level")
    elseif relax_type == "FFCF" then
        inst.braid:set_n_relax(-1, 1)
        inst.braid:set_n_relax(0, 0)
        print("    MGRIT uses F-FCF-Relaxation")
    elseif relax_type == "FCFF" then
        inst.braid:set_n_relax(-1, 0)
        inst.braid:set_n_relax(0, 1)
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
        end
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

    fine_level_inst = {}
    fine_level_inst.conv_check = util.xbraid.get_conv_check(desc,inst,desc.fine,fine_level_inst)
    fine_level_inst.solver = util.xbraid.get_solver(desc,inst,desc.fine,fine_level_inst)
    fine_level_inst.integrator = util.xbraid.get_integrator(desc,inst,desc.fine,fine_level_inst)
    driver:set_level_integrator(0,fine_level_inst.integrator)

    for level = 1, desc.max_level - 1 do
        coarse_level_inst = {}
        coarse_level_inst.conv_check = util.xbraid.get_conv_check(desc,inst,desc.coarse,coarse_level_inst)
        coarse_level_inst.solver = util.xbraid.get_solver(desc,inst,desc.coarse,coarse_level_inst)
        coarse_level_inst.integrator = util.xbraid.get_integrator(desc,inst,desc.coarse,coarse_level_inst)
        driver:set_level_integrator(level,coarse_level_inst.integrator)
        print("added to level lvl = " .. level)
    end
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


    driver:set_start_time(desc.time_interval.start_time)
    driver:set_end_time(desc.time_interval.end_time)
    driver:set_number_of_timesteps(desc.time.numtime)
    driver:set_time_values(desc.time.t_0, desc.time.t_end, desc.time.numtime)
    driver:set_max_levels(desc.max_level)
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





function util.xbraid.create_instance(desc,inst)
    braid = BraidExecutor(inst.communicator, inst.driver)
    inst.braid = braid
    util.xbraid.set_relax_and_cycle(desc, inst)


    if inst.prepared_time_hierarchy == nil then
        util.xbraid.create_time_hierarchy(desc,inst)
    end

    util.xbraid.create_hierarchy(desc,inst)
    print()
    return braid
end