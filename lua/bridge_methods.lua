util = {}
util.xbraid = {}

function util.xbraid.create_bdf_integrator(desc,inst)
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
    method = BDF_Integrator()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_order(desc.integrator.order)
    method:set_reassemble_threshold(desc.integrator.reassemble_threshold)
    return method
end


function util.xbraid.create_bdf_integrator_nl(desc,inst)
    if inst.domain_disc == nil then
        print("[ ERROR ] domain_disc must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator")
        exit()
    end
    if inst.nonlinear_solver == nil then
        print("[ ERROR ] nonlinear_solver must be set into *inst*")
        print("util.xbraid.create_theta_const_step_integrator")
        exit()
    end
    method = BDF_IntegratorNL()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_order(desc.integrator.order)
    return method
end

function util.xbraid.create_theta_const_step_integrator(desc,inst)
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
    method = ThetaConstStepIntegrator()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_num_steps(desc.integrator.num_steps)
    method:set_reassemble_threshold(desc.integrator.threshold)
    return method
end


function util.xbraid.create_theta_integrator_nl(desc,inst)
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
    method = ThetaIntegratorNL()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_order(desc.integrator.order)
    return method
end

function util.xbraid.create_theta_const_step_integrator_nl(desc,inst)
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
    method = ThetaConstStepIntegratorNL()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.nonlinear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_num_steps(desc.integrator.num_steps)
    return method
end

function util.xbraid.create_theta_single_timestep(desc,inst)
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
    method = ThetaSingleTimeStep()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_reassemble_threshold(desc.integrator.reassemble_threshold)
    return method
end
--[[
function util.xbraid.create_experimental_timestep(desc,inst)
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
    method = ExperimentalTimeStep()
    method:set_domain(inst.domain_disc)
    method:set_solver(inst.linear_solver)
    method:set_theta(desc.integrator.theta)
    method:set_rhs()
    method:set_jacobian()
    method:set_reassemble_threshold(desc.integrator.reassemble_threshold)
end
]]--
