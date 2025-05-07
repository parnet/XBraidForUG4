function util.xbraid.create_linear_conv_check(desc,inst)
    conv_check = ConvCheck()
    conv_check:set_maximum_steps(desc.linear_conv_check.maximum_steps)
    conv_check:set_reduction(desc.linear_conv_check.reduction)
    conv_check:set_minimum_defect(desc.linear_conv_check.minimum_defect)
    conv_check:set_verbose(desc.linear_conv_check.verbose)
    conv_check:set_supress_unsuccessful(desc.linear_conv_check.suppress_unsuccessful)
    return conv_check
end



function util.xbraid.create_gmg(desc,inst)
    local gmg = GeometricMultiGrid(approxSpace)
    gmg:set_discretization(domainDiscT)
    gmg:set_base_level(0)  -- was 1 in Cincy
    gmg:set_base_solver(superLU)  -- was baseLU in Cincy
    gmg:set_presmoother(preSmoother) --(jac)
    gmg:set_postsmoother(postSmoother)
    gmg:set_cycle_type("W") -- 1:V, 2:W -- "F"
    gmg:set_num_presmooth(2)
    gmg:set_num_postsmooth(2)
    gmg:set_rap(true)  -- mandatory, if set_stationary
    gmg:set_transfer(transfer)
    return gmg
end

function util.xbraid.create_bicgstab_gmg(desc,inst)
    method = BiCGStab()
    method:set_preconditioner(gmg)
    method:set_convergence_check(desc.linear_conv_check)
    return method
end

function util.xbraid.create_ls_gmg(desc,inst)
    method = LinearSolver()
    method:set_preconditioner(gmg)
    method:set_convergence_check(desc.linear_conv_check)
    return method
end

function util.xbraid.create_superlu(desc,inst)
    method =  LinearSolver()
    method:set_preconditioner(SuperLU())
    method:set_convergence_check(desc.linear_conv_check)
    return method
end

function util.xbraid.create_lu(desc,inst)
    method =  LinearSolver()
    method:set_preconditioner(LU())
    method:set_convergence_check(desc.linear_conv_check)
    return method
end

function util.xbraid.create_lu(desc,inst)

end

function util.xbraid.create_lu(desc,inst)

end

function util.xbraid.create_linear_solver(desc,inst)

end


function util.xbraid.create_newton_conv_check(desc,inst)
    local conv_check = ConvCheck()
    conv_check:set_maximum_steps(desc.newton_conv_check.maximum_steps) -- linear case!
    conv_check:set_reduction(desc.newton_conv_check.reduction)
    conv_check:set_minimum_defect(desc.newton_conv_check.minimum_defect) -- 1e-14
    conv_check:set_verbose(desc.newton_conv_check.verbose)
    conv_check:set_supress_unsuccessful(desc.newton_conv_check.suppress_unsuccessful)
    return conv_check
end


function util.xbraid.create_newton_solver(desc,inst)
    method = NewtonSolver()
    method:set_linear_solver(inst.linear_solver)
    method:set_convergence_check(inst.newton_conv_check)
    return method
end