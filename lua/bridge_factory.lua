
util = {}
util.xbraid = {}

function util.xbraid.create_limex_factory(desc,inst)
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

    method = LimexFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_error_estimator(inst.error_estimator)

    method:set_dt_max(desc.integrator.dt_max)
    method:set_dt_min(desc.integrator.dt_min)
    method:set_tol(desc.integrator.tol)
    method:set_level_factor(desc.integrator.level_factor)

    return method
end

function util.xbraid.create_linear_time_integrator_factory(desc,inst)
    method = LinearTimeIntegratorFactory()
    method:set_time_disc(inst.time_disc)
    method:set_solver(inst.linear_solver)
    return method
end


function util.xbraid.create_const_step_linear_time_integrator_factory(desc,inst)
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
    method = ConstStepLinearTimeIntegratorFactory()
    method:set_time_disc(inst.time_disc) -- inst.time_disc=ThetaTimeStep(domainDiscT)
    method:set_solver(inst.linear_solver)
    method:set_num_steps(desc.integrator.num_steps)
    return method
end

function util.xbraid.create_theta_integrator_factory(desc,inst)
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_integrator_factory")
        exit()
    end
    if inst.linear_solver == nil then
        print("[ ERROR ] linear_solver must be set into *inst*")
        print("util.xbraid.create_fixed_step_theta_factory")
        exit()
    end
    method = ThetaIntegratorFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.default_theta)
    for i, v in pairs(desc.integrator.theta) do
        method:set_level_theta(i,v)
    end
    return method
end

function util.xbraid.create_fixed_step_theta_factory(desc,inst)
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
    method = FixedStepThetaIntegratorFactory()
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



function util.xbraid.create_bdf_integrator_factory(desc,inst)
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
    method = BDF_IntegratorFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_order(desc.integrator.default_order)
    for i, v in pairs(desc.integrator.order) do
        method:set_level_order(i,v)
    end
    return method
end



function util.xbraid.create_simple_integrator_factory(desc,inst)
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
    method = SimpleIntegratorFactory()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_dt_min(desc.integrator.dt_min)
    method:set_dt_max(desc.integrator.dt_min)
    return method
end