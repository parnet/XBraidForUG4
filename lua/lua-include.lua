util = util or {}
util.xbraid = util.xbraid or {}







function message(meddelande)
    print(meddelande)
end







function dprint(meddelande)
    print(meddelande)
end






function typeprint(name, value)
    print(name)
    print("type = "..type(value))
    print(value)
    print("---------------------------")
end








-- =====================================================================================================================
--                                                   General
-- =====================================================================================================================
-- todo kollar om sektionen är klart !
-- todo test sektionen

function util.xbraid.split_communicator(desc,inst)
    message("<lua:util.xbraid.split_communicator>")
    local num_spatial_procs = util.GetParamNumber("--npx", 1,
            "number of spatial processors (must divide number of total_processors number)")
    local threads  = util.GetParamNumber("--threads", 1, "number threads for openmp")

    local num_world_ranks = NumProcs()
    local  space_time_communicator = SpaceTimeCommunicator()

    if num_world_ranks % num_spatial_procs == 0 then
        space_time_communicator:split(num_spatial_procs)
        local num_temporal_procs = num_world_ranks / num_spatial_procs;
        local num_spatial_procs = num_spatial_procs
        print("<xbraid:processor_configuration num_temporal_procs='"..num_temporal_procs.."' num_spatial_procs='" .. num_spatial_procs .. "' num_total_procs='" .. num_world_ranks.."' />")

    else
        space_time_communicator:split(1)
        print("[WARNING]    temporal processor does not divide total number of processors")
        print("<xbraid:processor_configuration num_temporal_procs='"..num_temporal_procs.."' num_spatial_procs='" .. num_spatial_procs .. "' num_total_procs='" .. num_world_ranks.."' />")

    end
    -- space_time_communicator:set_openmp(0,threads) -- todo openmp threads
    -- print("<xbraid:openmp threads='"..threads.."' />
    inst.communicator = space_time_communicator
    message("</lua:util.xbraid.split_communicator>")
    return space_time_communicator
end







function util.xbraid.cfactor_from_string(str_cfactor)
    message("<util.xbraid.cfactor_from_string>")
    local cfactor = {}
    local i = 1
    for s in str_cfactor:gmatch("[^_]+") do
        local factor = tonumber(s)
        print("level_"..i.." cfactor="..factor.." i="..i)
        cfactor[i] = factor
        i = i + 1

    end
    print("number of level by factor l="..i-1)
    message("</util.xbraid.cfactor_from_string>")
    return cfactor, i
end




function util.xbraid.nrelax_from_string(str_nrelax)
    message("<util.xbraid.nrelax_from_string>")
    local nrelax = {}
    local i = 1
    for s in str_nrelax:gmatch("[^_]+") do
        local factor = tonumber(s)
        print("level_"..i.." nrelax="..factor.." i="..i)
        nrelax[i] = factor
        i = i + 1

    end
    print("number of level by factor l="..i-1) -- todo kolla line
    message("</util.xbraid.nrelax_from_string>")
    return cfactor, i
end







function util.xbraid.numref_from_string(str_num_ref)
    message("<util.xbraid.numref_from_string>")
    local num_refs = {}
    local i = 1
    for s in str_num_ref:gmatch("[^_]+") do
        num_refs[i] = tonumber(s)
        i = i + 1
    end
    message("</util.xbraid.numref_from_string>")
    return num_refs, i
end






function util.xbraid.set_domain_disc(desc, inst, domain_disc)
    if domain_disc == nil then
        print("[ ERROR ]    domain_disc must be set in *parameter:domain_disc*")
        print("util.xbraid.set_domain_disc")
        exit()
    end
    message("<util.xbraid.set_domain_disc>")
    inst.domain_disc = domain_disc
    message("</util.xbraid.set_domain_disc>")
end






function util.xbraid.set_approx_space(desc,inst,approx_space)
    message("<util.xbraid.set_approx_space>")
    if approx_space == nil then
        print("[ ERROR ] approx_space are not allowed to be nil")
        print("util.xbraid.set_approx_space")
        exit()
    end
    inst.approx_space = approx_space
    message("</util.xbraid.set_approx_space>")
end







function util.xbraid.set_time_disc(desc,inst,time_disc)
    message("<util.xbraid.set_time_disc>")
    if time_disc == nil then
        print("[ ERROR ] time_disc are not allowed to be nil")
        print("util.xbraid.time_disc")
        exit()
    end
    inst.time_disc = time_disc
    message("</util.xbraid.set_time_disc>")
end







-- =====================================================================================================================
--                                                   Other
-- =====================================================================================================================
-- todo kollar om sektionen är klart !
-- todo test sektionen


function util.xbraid.create_spatial_grid_transfer(desc,inst)
    message("<util.xbraid.create_spatial_grid_transfer>")
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_spatial_grid_transfer")
        exit()
    end
    if inst.approx_space == nil then
        print("[ ERROR ] approx_space must be set into *inst*")
        print("util.xbraid.create_spatial_grid_transfer")
        exit()
    end
    local sgt = SpatialGridTransfer()
    sgt:set_transfer()
    sgt:set_domain(inst.domain_disc)
    sgt:set_approx_space(inst.approx_space)
    if inst.transfer ~= nil then
        sgt:set_transfer(inst.transfer)          -- todo userdata, string construction?
    elseif inst.transfer.prolongation ~= nil and inst.transfer.restriction ~= nil then
        sgt:set_prolongation(inst.prolongation)  -- todo userdata, string construction?
        sgt:set_restriction(inst.restriction)    -- todo userdata, string construction?
    else
        print("[ ERROR ] transfer or (prolongation and restriction) must be set into *inst*")
        print("util.xbraid.create_spatial_grid_transfer")
        exit()
    end
    message("</util.xbraid.create_spatial_grid_transfer>")
    return sgt
end







function util.xbraid.create_std_transfer(desc,inst)
    message("<util.xbraid.create_std_transfer>")
    local transfer = StdTransfer()
    transfer:enable_p1_lagrange_optimization(true)
    message("</util.xbraid.create_std_transfer>")
    return transfer
end







function util.xbraid.create_euclidian_norm(desc,inst)
    message("<util.xbraid.create_euclidian_norm>")
    local method = BraidEuclidianNorm();
    message("</util.xbraid.create_euclidian_norm>")
    return method
end







function util.xbraid.create_norm(desc,inst)
    message("<util.xbraid.create_norm>")
    if type(desc.norm) == "userdata" then
        return desc.norm

    elseif type(desc.norm) == "string" then
        if desc.norm == "l2" then
            local method = util.xbraid.create_euclidian_norm(desc,inst)
            message("</util.xbraid.create_norm>")
            return method
            -- todo inst.driver:set_norm_provider(method)
        else
            print("[ ERROR ]    Norm must set - given '", desc.norm, "' not implemented ")
            print("util.xbraid.create_norm")
            exit()
        end
    else
        print("[ ERROR ]    Norm must set")
        print("util.xbraid.create_norm")
        exit()
    end
    message("</util.xbraid.create_norm>")
end







function util.xbraid.create_start_value_initializer(desc,inst)
    local method = GridFunctionInitializer()
    return method
end






function util.xbraid.create_zero_initializer(desc,inst)
    local method = ZeroInitializer()
    return method
end







function util.xbraid.create_normal_initializer(desc,inst)
    local method = RandomValueInitializer()
    method:set_parameter_normal(desc.initializer.mean, desc.initializer.std)
    return method
end






function util.xbraid.create_uniform_initializer(desc,inst)
    local method = RandomValueInitializer()
    method:set_parameter_uniform(desc.initializer.min, desc.initializer.max)
    return method
end






function util.xbraid.create_initializer(desc,inst)
    if type(desc.initializer) == "userdata" then
        inst.initializer = desc.initializer
    elseif type(desc.initializer) == "table" then
        local name = desc.initializer.name
        if name == "default" then
            inst.initializer = util.xbraid.create_start_value_initializer(desc, inst)
        elseif name == "zero" then
            inst.initializer = util.xbraid.create_zero_initializer(desc, inst)
        elseif name == "uniform" then
            inst.initializer = util.xbraid.create_uniform_initializer(desc, inst)
        elseif name == "normal" then
            inst.initializer = util.xbraid.create_normal_initializer(desc, inst)
        else
            print("[WARNING] Invalid option for Initializer, using default type")
        end
    end
end







function util.xbraid.set_communicator(desc, inst, spc)
    if spc == nil then
        print("[ ERROR ]    spc must be set in *parameter:spc*")
        print("util.xbraid.set_communicator")
        exit()
    end
    message("<util.xbraid.set_communicator>")
    inst.communicator = spc
    message("</util.xbraid.set_communicator>")
end













-- =====================================================================================================================
--                                                   Integrator
-- =====================================================================================================================
-- todo kollar om sektionen är klart !
-- todo test sektionen


function util.xbraid.create_bdf_integrator(g_desc,g_inst,t_desc,t_inst)
    if g_inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator")
        exit()
    end
    if t_inst.linear_solver == nil then
        if g_inst.linear_solver == nil then
            print("[ ERROR ] linear_solver must be set into *inst*")
            print("util.xbraid.create_theta_const_step_integrator")
            exit()
        end
        t_inst.linear_solver = g_inst.linear_solver
    end

    local method = BDF_Integrator()

    method:set_domain(g_inst.domain_disc)

    method:set_solver(t_inst.linear_solver)
    method:set_order(t_desc.order)
    method:set_reassemble_threshold(t_desc.reassemble_threshold)

    return method
end







function util.xbraid.create_bdf_integrator_nl(g_desc,g_inst,t_desc,t_inst)
    if g_inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator")
        exit()
    end
    if t_inst == nil then
        if g_inst.nonlinear_solver == nil then
            print("[ ERROR ] nonlinear_solver must be set into *inst*")
            print("util.xbraid.create_theta_const_step_integrator")
            exit()
        end
        t_inst.nonlinear_solver = g_inst.nonlinear_solver
    end

    local method = BDF_IntegratorNL()

    method:set_domain(g_inst.domain_disc)

    method:set_solver(t_inst.nonlinear_solver)
    method:set_order(t_desc.order)
    return method
end







function util.xbraid.create_theta_const_step_integrator(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator")
        exit()
    end
    local method = ThetaConstStepIntegrator()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_num_steps(desc.integrator.num_steps)
    method:set_reassemble_threshold(desc.integrator.threshold)
    return method
end







function util.xbraid.create_theta_integrator_nl(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    print(t_desc)
    message("<util.xbraid.create_theta_integrator_nl>")
    if g_inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator_nl")
        exit()
    end
    print(t_desc)
    if t_desc.nonlinear_solver == nil then -- todo more search options g_desc, g_inst ?
        print("[ ERROR ] nonlinear_solver must be set into *t_desc*")
        print("util.xbraid.create_theta_const_step_integrator_nl")
        exit()
    end

    local method = ThetaIntegratorNL()
    method:set_domain(g_inst.domain_disc)
    method:set_solver(t_desc.nonlinear_solver)

    method:set_theta(t_desc.theta)
    message("</util.xbraid.create_theta_integrator_nl>")
    return method
end






function util.xbraid.create_theta_const_step_integrator_nl(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator_nl")
        exit()
    end
    if inst.nonlinear_solver == nil then
        print("[ ERROR ] nonlinear_solver must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator_nl")
        exit()
    end
    local method = ThetaConstStepIntegratorNL()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_num_steps(desc.integrator.num_steps)
    return method
end






function util.xbraid.create_theta_single_timestep(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_single_timestep")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_theta_single_timestep")
        exit()
    end
    local method = ThetaSingleTimeStep()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_reassemble_threshold(desc.integrator.reassemble_threshold)
    return method
end






--[[
function util.xbraid.create_experimental_timestep(desc,inst)
-- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_experimental_timestep")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_experimental_timestep")
        exit()
    end
    local method = ExperimentalTimeStep()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_rhs()
    method:set_jacobian()
    method:set_reassemble_threshold(desc.integrator.reassemble_threshold)
end
]]--


function util.xbraid.create_limex_factory(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.nonlinear_solver == nil then
        print("[ ERROR ] nonlinear_solver must be set into *inst*")
        print("util.xbraid.create_limex_factory")
        exit()
    end
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_limex_factory")
        exit()
    end
    if inst.error_estimator == nil then
        print("[ ERROR ] error_estimator must be set into *inst*")
        print("util.xbraid.create_limex_factory")
        exit()
    end

    local method = LimexFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_error_estimator(inst.error_estimator)

    method:set_dt_max(desc.integrator.dt_max)
    method:set_dt_min(desc.integrator.dt_min)
    method:set_tol(desc.integrator.tol)
    method:set_level_factor(desc.integrator.level_factor)

    return method
end







function util.xbraid.create_linear_time_integrator_factory(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    local method = LinearTimeIntegratorFactory()
    method:set_time_disc(inst.time_disc)
    method:set_solver(inst.linear_solver)
    return method
end







function util.xbraid.create_const_step_linear_time_integrator_factory(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.time_disc == nil then
        print("[ ERROR ] time_disc must be set into *inst*")
        print("util.xbraid.create_const_step_linear_time_integrator_factory")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_const_step_linear_time_integrator_factory")
        exit()
    end
    local method = ConstStepLinearTimeIntegratorFactory()
    method:set_time_disc(inst.time_disc) -- inst.time_disc=ThetaTimeStep(domainDiscT)
    method:set_solver(inst.linear_solver)
    method:set_num_steps(desc.integrator.num_steps)
    return method
end






function util.xbraid.create_theta_integrator_factory(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_integrator_factory")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_theta_integrator_factory")
        exit()
    end
    local method = ThetaIntegratorFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.default_theta)
    for i, v in pairs(desc.integrator.theta) do
        method:set_level_theta(i,v)
    end
    return method
end






function util.xbraid.create_fixed_step_theta_factory(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_fixed_step_theta_factory")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_fixed_step_theta_factory")
        exit()
    end
    local method = FixedStepThetaIntegratorFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.default_theta)
    for i, v in pairs(desc.integrator.order) do
        method:set_level_theta(i,v)
    end
    for i, v in pairs(desc.integrator.num_steps) do
        method:set_level_num_steps(i,v)
    end
    method:create_time_integrator()
    method:create_level_time_integrator()
    return method
end






function util.xbraid.create_bdf_integrator_factory(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_bdf_integrator_factory")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_bdf_integrator_factory")
        exit()
    end
    local method = BDF_IntegratorFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_order(desc.integrator.default_order)
    for i, v in pairs(desc.integrator.order) do
        method:set_level_order(i,v)
    end
    return method
end






function util.xbraid.create_simple_integrator_factory(g_desc,g_inst,t_desc,t_inst)
    -- todo change desc -> g_desc
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_simple_integrator_factory")
        exit()
    end
    if inst.nonlinear_solver == nil then
        print("[ ERROR ] nonlinear_solver must be set into *inst*")
        print("util.xbraid.create_simple_integrator_factory")
        exit()
    end
    local method = SimpleIntegratorFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_dt_min(desc.integrator.dt_min)
    method:set_dt_max(desc.integrator.dt_min)
    -- todo or delete method:set_reduction_factor(desc.integrator.reduction_factor) -- default=0.2
    return method
end







function util.xbraid.create_target_integrator(g_desc,g_inst, t_desc, t_inst)
    message("<util.xbraid.create_target_integrator>")
    local method
    if t_desc.name == "ThetaSingleStep" then
        method = util.xbraid.create_theta_single_timestep(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "ThetaIntegrator" then
        print("not implemented - ThetaIntegrator - use ThetaSingleStep instead")
        print("util.xbraid.create_target_integrator")
        exit()

    elseif t_desc.name == "ThetaIntegratorNL" then
        method = util.xbraid.create_theta_integrator_nl(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "ThetaConstStepIntegrator" then
        method = util.xbraid.create_theta_const_step_integrator(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "ThetaConstStepIntegratorNL" then
        method = util.xbraid.create_theta_const_step_integrator_nl(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "ExperimentalTimeStep" then
        print("not implemented yet - ExperimentalTimeStep")
        print("util.xbraid.create_target_integrator")
        exit()

    elseif t_desc.name == "BDF_Integrator" then
        method = util.xbraid.create_bdf_integrator(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "BDF_IntegratorNL" then
        method = util.xbraid.create_bdf_integrator_nl(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "ThetaIntegratorFactory" then
        g_desc.factory = true
        method = util.xbraid.create_theta_integrator_factory(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "SimpleIntegratorFactory" then
        g_desc.factory = true
        method = util.xbraid.create_simple_integrator_factory(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "LinearTimeIntegratorFactory" then
        g_desc.factory = true
        method = util.xbraid.create_linear_time_integrator_factory(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "LimexFactory" then
        g_desc.factory = true
        method = util.xbraid.create_limex_factory(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "FixedStepThetaIntegratorFactory" then
        g_desc.factory = true
        method = util.xbraid.create_fixed_step_theta_factory(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "ConstStepLinearTimeIntegratorFactory" then
        g_desc.factory = true
        method = util.xbraid.create_const_step_linear_time_integrator_factory(g_desc,g_inst,t_desc,t_inst)

    elseif t_desc.name == "BDF_IntegratorFactory" then
        g_desc.factory = true
        method = util.xbraid.create_bdf_integrator_factory(g_desc,g_inst,t_desc,t_inst)

    else
        print("[ ERROR ]   Name not specified")
        print("util.xbraid.create_target_integrator")
        exit()
    end
    message("</util.xbraid.create_target_integrator>")
    t_inst.integrator = method
    return method
end




function util.xbraid.create_integrator(desc, inst)
    message("<util.xbraid.create_integrator>")
    message("<default_integrator>")
    inst.integrator = {}
    util.xbraid.create_target_integrator(desc, inst, desc.integrator, inst.integrator)
    message("</default_integrator>")
    if desc.level ~= nil then
        inst.level = {}
        for level, k in pairs(desc.level) do
            inst.level[level] = {}
            print("<util.xbraid.create_integrator level=\""..level.."\">")
            util.xbraid.create_level_integrator(desc, inst,desc.level[level], inst.level[level])
            print("</util.xbraid.create_integrator level=\""..level.."\">")
        end
    end
    message("</util.xbraid.create_integrator>")
    return method
end





-- =====================================================================================================================
--                                                   Driver
-- =====================================================================================================================
-- todo kollar om sektionen är klart !
-- todo test sektionen

function util.xbraid.init_driver_base(desc,inst,method)
    message("<util.xbraid.init_driver_base>")
    if inst.spatial_norm == nil then
        print("[ ERROR ]    spatial_norm must be set in *inst*")
        print("util.xbraid.init_driver_base")
        exit()
    end
    if inst.domain_disc == nil then
        print("[ ERROR ]    spatial_norm must be set in *inst*")
        print("util.xbraid.init_driver_base")
        exit()
    end
    -- method:set_start_time()
    -- method:set_end_time()
    -- method:set_number_of_timesteps()
    method:set_time_values(desc.time_interval.start_time,
            desc.time_interval.end_time,
            desc.time_interval.time_steps)
    -- method:set_start_vector() is set by executor
    method:set_norm_provider(inst.spatial_norm)
    if inst.process_observer ~= nil then
        method:attach_xbraid_observer(inst.process_observer)
    end
    if inst.observer ~= nil then
        method:attach_observer(inst.observer)
    end
    method:set_max_levels(desc.hierarchy.max_levels)
    method:set_domain(inst.domain_disc)

    if inst.initializer ~= nil then
        method:set_initializer(inst.initializer) -- optional default StartValueInitializer
    end
    -- method:init() -- todo kollar när det måste användas, log måste sättas först
    message("</util.xbraid.init_driver_base>")
end







function util.xbraid.create_basic_driver(desc, inst)
    message("<util.xbraid.create_basic_driver>")
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set in *inst*")
        print("util.xbraid.create_basic_driver")
        exit()
    end
    if inst.integrator == nil then
        print("[ ERROR ] default_integrator must be set in *inst*")
        print("util.xbraid.create_basic_driver")
        exit()
    end
    local method = BasicDriver()
    util.xbraid.init_driver_base(desc,inst,method)

    method:set_domain(inst.domain_disc)
    method:set_default_integrator(inst.integrator.integrator)
    if desc.level_config == "fine_coarse" then
        method:set_level_integrator() --  todo
    end
    for level in desc.max_levels do
        dprint("[ TODO  ] ".. level)
        method:set_level_integrator() --  todo
    end


    if inst.transfer ~= nil then
        method:set_spatial_grid_transfer(inst.transfer)
        method:set_level_num_ref() -- todo
    end
    message("</util.xbraid.create_basic_driver>")
    return method
end








function util.xbraid.create_braid_integrator(desc, inst)
    message("<util.xbraid.create_braid_integrator>")
    local method = BraidIntegrator()
    util.xbraid.init_driver_base(desc,inst,method)


    method:set_ref_factor(desc.ref_factor)
    method:set_threshold(desc.threshold)

    method:set_default_integrator(inst.integrator.integrator)
    for level in levels do
        method:set_integrator(level, inst.integrator)
    end
    message("</util.xbraid.create_braid_integrator>")
    return method
end







function util.xbraid.create_nl_integrator(desc, inst)
    message("<util.xbraid.create_nl_integrator>")
    local method = BraidNLIntegrator()
    util.xbraid.init_driver_base(desc,inst,method)

    method:set_default_integrator(inst.integrator.integrator)
    -- todo for level in levels do
    -- todo     method:set_integrator(level, inst.level_level_integrator[level])
    -- todo  end
    -- todo (empty c++ method) method:set_conv_check(inst.conv_check)
    -- todo (unused c++ variable) method:set_threshold(desc.threshold)
    -- todo (unused c++ variable/method) method:set_tol(desc.tol)
    message("</util.xbraid.create_nl_integrator>")
    return method
end








function util.xbraid.create_integrator_factory(desc, inst)
    message("<util.xbraid.create_integrator_factory>")
    local method = BraidIntegratorFactory()

    util.xbraid.init_driver_base(desc, inst, method)

    method:set_default_integrator(inst.integrator.integrator)
    for level in levels do
        method:set_coarse_time_integrator(level, inst.coarse_integrator)
    end
    message("</util.xbraid.create_integrator_factory>")
    return method
end









function util.xbraid.create_driver(desc, inst)
    message("<util.xbraid.create_driver>")
    if desc.driver == nil then
        if inst.driver == nil then
            print("[ ERROR ]    driver has to be set in *desc* or *inst*")
            print("util.xbraid.create_driver")
            exit()
        end
        print("[ INFO  ]    using inst.driver instance - empty desc")
        message("</util.xbraid.create_driver>")
        return inst.driver
    end
    local method

    if type(desc.driver) == "userdata" then
        message("[ INFO  ]    using desc.driver instance (userdata)")
        message("</util.xbraid.create_driver>")
        inst.driver = desc.driver
        return desc.driver
    else -- type == assuming table
        if desc.driver.name == "BasicDriver" then
            method = util.xbraid.create_basic_driver(desc,inst)
            message("</util.xbraid.create_driver>")

        elseif desc.driver.name == "Integrator" then
            method = util.xbraid.create_braid_integrator(desc,inst)
            message("</util.xbraid.create_driver>")

        elseif desc.driver.name == "Nonlinear" then
            method = util.xbraid.create_nl_integrator(desc,inst)
            message("</util.xbraid.create_driver>")

        elseif desc.driver.name == "Factory" then
            method = util.xbraid.create_integrator_factory(desc,inst)
            message("</util.xbraid.create_driver>")

        else
            print("[ ERROR ]    ".. desc.driver.name .." driver name is unknwon")
            print("util.xbraid.create_driver")
            exit()
        end
    end
    message("</util.xbraid.create_driver>")
    inst.driver = method
    return method
end







-- =====================================================================================================================
--                                                    Observer
-- =====================================================================================================================
-- todo kollar om sektionen är klart !
-- todo test sektionen

function util.xbraid.create_time_integrator_observer_collector(desc,inst)
    local method = TimeIntegratorObserverCollector()
    method:attach_observer()
    return method
end







function util.xbraid.create_xBraid_time_integrator_observer_collector(desc,inst)
    local method = XBraidTimeIntegratorObserverCollector()
    method:attach_observer()
    method:attach_common_observer()
    return method
end







function util.xbraid.create_vtk_observer(observer_desc,inst)
    if inst.domain_disc == nil then
        print("[ ERROR ] vtk must be set into *inst*")
        print("util.xbraid.create_vtk_process_observer")
        exit()
    end
    if desc.filename == nil then
        print("[WARNING] filename should be set into *desc*")
        print("util.xbraid.create_vtk_process_observer")
        desc.filename = "vtk_file"
    end
    local method = VTK_Observer(inst.vtk, desc.filename)
    return method
end






function util.xbraid.create_matlab_observer(desc,inst)
    local method = MATLAB_Observer()
    return method
end






function util.xbraid.create_vtk_process_observer(desc,inst)
    if inst.domain_disc == nil then
        print("[ ERROR ] vtk must be set into *inst*")
        print("util.xbraid.create_vtk_process_observer")
        exit()
    end
    if desc.filename == nil then
        print("[WARNING] filename should be set into *desc*")
        print("util.xbraid.create_vtk_process_observer")
        desc.filename = "vtk_file"
    end
    local method = VTK_ProcessObserver(inst.vtk, desc.filename)
    return method
end






function util.xbraid.create_eval_observer(desc,inst)
    local method = EvalObserver()
    if inst.domain_disc == nil then
        print("[ ERROR ] domain disc must be set into *inst*")
        print("util.xbraid.create_eval_observer")
        exit()
    end
    method:set_domain(inst.domain_disc)
    method:set_filename(desc.filename)
    method:set_generator_component(desc.component) -- const char *
    method:set_vector_generator(desc.generator) -- smartpointer to userdata
    method:set_relative(desc.relative) -- bool
    return method
end





function util.xbraid.create_single_observer(desc,inst)
    local observer
    if desc.name == "EvalObserver" then
        observer = util.xbraid.create_eval_observer(desc,inst)
    elseif desc.name == "VTKProcessObserver" then
        observer = util.xbraid.create_vtk_process_observer(desc,inst)
    elseif desc.name == "MATLABObserver" then
        observer = util.xbraid.create_matlab_observer(desc,inst)
    elseif desc.name == "VTKObserver" then
        observer = util.xbraid.create_vtk_observer(desc,inst)
    else
        print("Not type named")
        print("util.xbraid.create_single_observer")
        exit()
    end
    return observer
end





-- todo split process observer and observer
function util.xbraid.create_observer(desc,inst)
    if type(desc.observer) == "userdata" then
        inst.observer = desc.observer
    elseif type(desc.observer) == "table" then
        if desc.observer.name == "ProcessObserverCollector" then
            local observer = util.xbraid.create_xBraid_time_integrator_observer_collector(desc,inst)
            inst.observer = observer
            for k, v in pairs(desc.observer.observers) do
                observer = util.xbraid.create_single_observer(v)
                observer.attach_common_observer(observer)
            end
            for k, v in pairs(desc.observer.process_observers) do
                observer = util.xbraid.create_single_observer(v)
                observer.attach_observer(observer)
            end
        elseif  desc.observer.name == "ObserverCollector" then
            local observer = util.xbraid.create_time_integrator_observer_collector(desc,inst)
            inst.observer = observer
            for k, v in pairs(desc.observer.observers) do
                observer = util.xbraid.create_single_observer(v)
                observer.attach_observer(observer)
            end
        else
            local observer = util.xbraid.create_single_observer(desc.observer.name)
            inst.observer = observer
        end
    end
end




function util.xbraid.create_logger(desc,inst)
    message("<util.xbraid.create_logger>")
    if inst.communicator == nil then
        print(" [ ERROR ]     communicator not set in *inst*")
        print("util.xbraid.create_logger")
        exit()
    end

    local log_job = Paralog()
    log_job:set_comm(inst.communicator)
    log_job:set_filename(desc.logger or "job") -- todo userdata or string
    log_job:init()
    inst.logger = log_job
    message("</util.xbraid.create_logger>")
    return logger
end



-- =====================================================================================================================
--                                                    Executor
-- =====================================================================================================================
-- todo kollar om sektionen är klart !
-- todo test sektionen




function util.xbraid.create_grid_hierarchy(desc,inst)
    message("<util.xbraid.create_grid_hierarchy>")
    -- level
    -- todo level    inst.xbraid:set_n_relax(level, desc.n_relax[level])
    -- todo level    inst.xbraid:set_c_factor(level, desc.c_factor[level])
    -- todo level    inst.xbraid:set_c_relax_weight(level, desc.c_relax_weight[level])
    inst.xbraid:set_c_relax_weight(-1, desc.cf_weight)
    -- method:set_max_levels(desc.max_levels)

    print("    Create time hierarchy!: ")
    if desc.default_cfactor ~= nil then
        inst.xbraid:set_c_factor(-1, desc.default_cfactor)
    end



    if type(desc.cfactor) == "string" then
        print(desc.cfactor)
        desc.cfactor, max_level_cfactor = util.xbraid.cfactor_from_string(desc.cfactor)
        print(desc.cfactor)
    end

    if type(desc.cfactor) == "table" then
        for key, value in pairs(desc.cfactor) do
            inst.xbraid:set_c_factor(key - 1, value)
        end
    end

    if type(desc.level_num_ref) == "string" then
        print(desc.level_num_ref)
        desc.level_num_ref = util.xbraid.numref_from_string(desc.level_num_ref)
        print(desc.level_num_ref)
    else
        print("desc.level_num_ref is not string")
    end

    if type(desc.level_num_ref) == "table" then
        for key, value in pairs(desc.level_num_ref) do
            -- todo process spatial coarsening
            print(" [WARNING] level num ref factor is discarded")
        end
    end

    desc_max_level = desc.max_level

    fnum_time = desc.time_interval.time_steps
    print("    \tlvl  | \tcfac | \tnum-time | \tnum-ref")
    print("    \t--------------------------------------")
    base_reached = false
    for i = 1, #desc.cfactor do
        cfac = desc.cfactor[i] -- todo or default if not set
        numref = desc.level_num_ref[i] -- todo or default if not set or last refinement level
        fnum_time_p = fnum_time / cfac  -- todo change to real coarsening factor for integer divion
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
    print()

    message("</util.xbraid.create_grid_hierarchy>")

end



function util.xbraid.create_braid_executor(desc,inst)
    message("<util.xbraid.create_braid_executor>")
    if inst.driver == nil then
        print("[ ERROR ]    driver must be set first in *inst*")
        print("util.xbraid.create_braid_executor")
        exit()
    end
    if inst.communicator == nil then
        print("[ ERROR ]    communicator must be set first in *inst*")
        print("util.xbraid.create_braid_executor")
        exit()
    end

    local method = BraidExecutor(inst.communicator, inst.driver)
    inst.xbraid = method
    if inst.initializer ~= nil then -- optional initializer default StartValueInitializer
        method:set_initializer(inst.initializer)
    end
    if inst.driver == nil then
        -- todo
        print("driver is not set")
    end


    if inst.logger == nil then
        -- todo
        print("logger is not set")
    end

    -- todo ø error - method:set_driver(inst.driver)
    method:set_parallel_logger(inst.logger) -- todo set app

    -- method:set_norm_provider() delegates to driver
    method:set_residual(desc.use_residual or false)
    util.xbraid.create_grid_hierarchy(desc,inst)

    -- todo ø print("__ min_coarse=".. desc.min_coarse)
    -- todo ø method:set_min_coarse(desc.min_coarse)

    print("__ desc.conv_check.max_iterations=".. desc.conv_check.max_iterations)
    method:set_max_iterations(desc.conv_check.max_iterations)
    if desc.conv_check.absolut_tol ~= 0 then
        print("__ desc.conv_check.absolut_tol=".. desc.conv_check.absolut_tol)
        method:set_absolute_tol(desc.conv_check.absolut_tol)
    end
    if desc.conv_check.relative_tol ~= 0 then
        print("__ desc.conv_check.relative_tol=".. desc.conv_check.relative_tol)
        method:set_relative_tol(desc.conv_check.relative_tol)
    end

    print("__ desc.temporal_norm=".. desc.temporal_norm)
    method:set_temporal_norm(desc.temporal_norm)

    if desc.sequential ~= nil then
        print("__ desc.sequential=".. desc.sequential)
        method:set_sequential(desc.sequential)
    end

    if desc.store_values ~= nil then
        print("__ desc.store_values=".. desc.store_values)
        method:set_store_values(desc.store_values)
    end

    if desc.spatial_coarsen_and_refine ~= nil then
        method:set_spatial_coarsen_and_refine(desc.spatial_coarsen_and_refine)
    end

    if desc.refine ~= nil then
        print("__ value=".. desc.refine)
        method:set_refine(desc.refine)

        print("__ desc.max_refinements=".. desc.max_refinements)
        method:set_max_refinements(desc.max_refinements)
    end

    print("__ desc.access_level=".. desc.access_level)
    method:set_access_level(desc.access_level)

    print("__ desc.print_level=".. desc.print_level)
    method:set_print_level(desc.print_level)

    print("__ desc.print_file=".. desc.print_file)
    method:set_print_file(desc.print_file)
    typeprint("desc.skip_downcycle_work",desc.skip_downcycle_work)
    if type(desc.skip_downcycle_work) == "number" then
        method:set_skip_downcycle_work(desc.skip_downcycle_work == 1)
    else
        method:set_skip_downcycle_work(desc.skip_downcycle_work)
    end

    local processed_relaxation = false
    if desc.relaxation ~= nil then
        if type(desc.relaxation) == "string" then
            if desc.relaxation == "F" then
                method:set_n_relax(-1, 0)
                processed_relaxation = true
            elseif desc.relaxation == "FCF" then
                method:set_n_relax(-1, 1)
                processed_relaxation = true
            elseif desc.relaxation == "FFCF" then
                method:set_n_relax(0, 0)
                method:set_n_relax(-1, 1)
                processed_relaxation = true
            end
        else
            message("invalid option for desc.relaxation")
            message("util.xbraid.create_braid_executor")
            exit()
        end
    end

    if not processed_relaxation then
        if desc.n_relax_default ~= nil then
            method:set_n_relax(-1, desc.n_relax_default)
            processed_relaxation = true
        end

        if type(desc.n_relax) == "string" then
            desc.n_relax = util.xbraid.nrelax_from_string(desc.n_relax)
        end

        if type(desc.n_relax) == "table" then
            for key, value in pairs(desc.n_relax) do
                method:set_n_relax(key - 1, value)
            end
            processed_relaxation = true
        end
        if not processed_relaxation then
            message("relaxation method not set")
            message("util.xbraid.create_braid_executor()")
            exit()
        end
    end

    if desc.n_relax then
    method:set_n_relax(level, number)
    end

    if desc.cycle.cycle_type == "F" then -- todo
    print("__ FMG=YES")
    method:set_cycle_fmg()

    print("__ desc.cycle.nfmg=".. desc.cycle.nfmg)
    method:set_cycle_nfmg(desc.cycle.nfmg)

    print("__ desc.cycle.nfmgv=".. desc.cycle.nfmgv)
    method:set_cycle_nfmgv(desc.cycle.nfmgv)
    print("MGRIT -> using F cycle")
    else
    print("MGRIT -> using V cycle")
    end

    --    print("__ desc.sync=".. desc.sync)
    if desc.sync then
    method:set_sync()
    end

    --print("__ desc.increase_max_levels=".. desc.increase_max_levels)
    if desc.increase_max_levels ~= nil then
    method:set_increase_max_levels()
    end

    --print("__ desc.relax_only_cg=".. desc.relax_only_cg)
    if desc.relax_only_cg ~= nil then
    method:set_relax_only_cg(desc.relax_only_cg)
    end

    if desc.agg_c_factor then
    print("__ desc.agg_c_factor=".. desc.agg_c_factor)
    method:set_agg_c_factor(desc.agg_c_factor)
    end

    if desc.periodic ~= nil then
    print("__ desc.periodic=".. desc.periodic)
    method:set_periodic(desc.periodic)
    end

    if desc.final_fc_relax ~= nil  then
    print("__ desc.final_fc_relax=".. desc.final_fc_relax)
    method:set_final_fc_relax(desc.final_fc_relax)
    end

    if desc.reverted_ranks ~= nil then
    print("__ desc.reverted_ranks=".. desc.reverted_ranks)
    method:set_reverted_ranks(desc.reverted_ranks)
    end

    if desc.richardson_estimation ~= nil then -- ø todo extrapolation and estimation are two different concepts that shares local_order
    print("__ desc.richardson_estimation=".. desc.richardson_estimation)
    print("__ desc.use_extrapolation=".. desc.use_extrapolation)
    print("__ desc.local_order=".. desc.local_order)
    method:set_richardson_estimation(desc.richardson_estimation, desc.use_extrapolation, desc.local_order)
    end

    if desc.file_io_level ~= nil then
    print("__ desc.file_io_level=".. desc.file_io_level)
    method:set_file_io_level(desc.file_io_level)
    end

    if desc.t_points_cutoff ~= nil then
    print("__ desc.t_points_cutoff=".. desc.t_points_cutoff)
    method:set_t_points_cutoff(desc.t_points_cutoff)
    end

    if desc.full_residual_norm then
    print("__ desc.full_residual_norm=YES")
    method:set_full_residual_norm()
    end

    -- todo ø print("__ value=".. inst.value)
    -- todo ø method:set_time_grid() -- todo empty method

    message("</util.xbraid.create_braid_executor>")
    return method
    end







    function util.xbraid.create_instance(desc,inst)
    message("<util.xbraid.create_instance>")
    if inst.xbraid == nil then
    util.xbraid.create_braid_executor(desc, inst)
    end
    message("</util.xbraid.create_instance>")
    return inst.xbraid
    end



    function util.xbraid.apply(desc, inst, u_start)
    print("<util.xbraid.apply>")
    local timer = BraidTimer()
    timer:start()
    inst.xbraid:apply(u_start, desc.time_interval.end_time, u_start, desc.time_interval.start_time)
    timer:stop()
    time = timer:get()
    print("<execution time=\""..time.."\" />")
    print("</util.xbraid.apply>")
    end



    -- =====================================================================================================================
    --
    -- =====================================================================================================================