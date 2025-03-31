-- ----------------------------------------------------------------------------

util = util or {}
util.xbraid = util.xbraid or {}
util.poro = util.poro or {}


-- ----------------------------------------------------------------------------
-- Conv Check
-- ----------------------------------------------------------------------------
function util.poro.create_conv_check(approxSpace, dim)
    local p0 = 1.0
    local cmpConvCheck = CompositeConvCheck(approxSpace)
    cmpConvCheck:set_component_check("ux", p0 * 1e-14, 1e-6)
    cmpConvCheck:set_component_check("uy", p0 * 1e-14, 1e-6)
    if (dim == 3) then -- todo dim from approx space?
        cmpConvCheck:set_component_check("uz", p0 * 1e-14, 1e-6)
    end
    cmpConvCheck:set_component_check("p", p0 * 1e-14, 1e-6)
    cmpConvCheck:set_maximum_steps(200) -- coarse = coarse ?
    cmpConvCheck:set_verbose(true)
    cmpConvCheck:set_supress_unsuccessful(true)
    return cmpConvCheck
end

function util.xbraid.create_conv_check()
    local convCheck = ConvCheck()
    convCheck:set_maximum_steps(200)
    convCheck:set_reduction(1e-8)
    convCheck:set_minimum_defect(1e-14)
    convCheck:set_verbose(true)
    convCheck:set_supress_unsuccessful(false)
    return convCheck
end


-- ----------------------------------------------------------------------------
-- Smoother
-- ----------------------------------------------------------------------------
function util.poro.create_uzawa(sSchurCmp, aiForward, aiSchur, aiBackward, uzawaSchurUpdateOp, uzawaSchurWeight)
    local uzawa = UzawaBase(sSchurCmp)
    local weight = uzawaSchurWeight or 1.0
    if (aiForward) then
        uzawa:set_forward_iter(aiForward)
    end
    if (aiSchur) then
        uzawa:set_schur_iter(aiSchur)
    end
    if (aiBackward) then
        uzawa:set_backward_iter(aiBackward)
    end
    uzawa:set_schur_operator_update(uzawaSchurUpdateOp, weight)
    return uzawa
end

function util.poro.create_smoother(smootherID)
    if smootherID == "uzawa1" then
        print("Using uzawa #1 as smoothing method")
        preSmoother = util.xbraid.create_uzawa("p", gs, Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = util.xbraid.create_uzawa("p", nil, Jacobi(0.66), bgs, uzawaSchurUpdateOp, uzawaWeight)

    elseif smootherID == "uzawa2" then
        print("Using uzawa #2 as smoothing method")
        preSmoother = util.xbraid.create_uzawa("p", gs, SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = util.xbraid.create_uzawa("p", nil, SymmetricGaussSeidel(), bgs, uzawaSchurUpdateOp, uzawaWeight)

    elseif smootherID == "uzawa3" then
        print("Using uzawa #3 as smoothing method")
        preSmoother = util.xbraid.create_uzawa("p", SymmetricGaussSeidel(), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = util.xbraid.create_uzawa("p", nil, Jacobi(0.66), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

    elseif smootherID == "uzawa4" then
        print("Using uzawa #4 as smoothing method")
        preSmoother = util.xbraid.create_uzawa("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = util.xbraid.create_uzawa("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

    elseif smootherID == "uzawa5" then
        print("Using uzawa #5 as smoothing method")
        preSmoother = util.xbraid.create_uzawa("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = util.xbraid.create_uzawa("p", nil, SymmetricGaussSeidel(), bgs, uzawaSchurUpdateOp, uzawaWeight)

    elseif smootherID == "uzawa6" then
        print("Using uzawa #6 as smoothing method")
        preSmoother = util.xbraid.create_uzawa("p", gs, SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
        postSmoother = util.xbraid.create_uzawa("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

    end
    return preSmoother, postSmoother
end


-- ----------------------------------------------------------------------------
-- Solver
-- ----------------------------------------------------------------------------

function util.xbraid.create_lsolver(selector, approxSpace, domainDiscT, preSmoother, postSmoother)
    local transfer = StdTransfer()
    transfer:enable_p1_lagrange_optimization(true)

    local superLU = SuperLU() --LU()
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

    convCheck = util.xbraid.create_conv_check()
    local solver = {}

    solver["GMG"] = LinearSolver()
    solver["GMG"]:set_preconditioner(gmg) -- gmg, dbgIter
    solver["GMG"]:set_convergence_check(convCheck)

    solver["SuperLU"] = LinearSolver()
    solver["SuperLU"]:set_preconditioner(SuperLU())
    solver["SuperLU"]:set_convergence_check(convCheck)

    solver["LU"] = LinearSolver()
    solver["LU"]:set_preconditioner(LU())
    solver["LU"]:set_convergence_check(convCheck)

    solver["GMGKrylov"] = BiCGStab()
    solver["GMGKrylov"]:set_preconditioner(gmg) -- gmg, dbgIter
    solver["GMGKrylov"]:set_convergence_check(convCheck) -- convCheck

    return solver[selector]
end

function util.xbraid.create_nlsolver(lsolver)
    local newtonCheck = ConvCheck()
    newtonCheck:set_maximum_steps(1) -- linear case!
    newtonCheck:set_reduction(9.999999999999e-1)  -- todo warum fast 1? 5e-6
    newtonCheck:set_minimum_defect(1e-14) -- 1e-14
    newtonCheck:set_verbose(true)
    newtonCheck:set_supress_unsuccessful(true)

    local newtonSolver = NewtonSolver()
    newtonSolver:set_linear_solver(lsolver)
    newtonSolver:set_convergence_check(newtonCheck)
    return newtonSolver
end
-- ----------------------------------------------------------------------------
-- TimeIntegrator
-- ----------------------------------------------------------------------------

function util.xbraid.create_NLThetaIntegrator()
    if num_steps == 1 then
        integrator = NLThetaIntegrator()
        integrator:set_domain(domain)
        integrator:set_solver(nlsolver)
        integrator:set_theta(theta)
        return integrator
    else
        integrator = NLFixedStepThetaIntegrator()
        integrator:set_domain(domain)
        integrator:set_solver(nlsolver)
        integrator:set_theta(theta)
        integrator:set_num_steps(num_steps)
        return integrator
    end
end


function util.xbraid.create_BDF_IntegratorFactory()
    integrator = BDF_IntegratorFactory()
    integrator:set_domain(domainDiscT)
    integrator:set_solver(lsolver)
    for i = 1, max_level do
        integrator:set_level_order(i-1, orderOrTheta)
    end
    return integrator
end

function util.xbraid.create_Limex_IntegratorFactory()
    integrator = LimexFactory()
    integrator:set_domain_disc(domainDiscT)
    integrator:set_solver(nlsolver)
    integrator:set_error_estimator(biotErrorEst)
    integrator:set_tol(1e-3)
    integrator:set_dt_min(1e-20)
    return integrator
end

function util.xbraid.create_LinearTimeIntegratorFactory()
    integrator = LinearTimeIntegratorFactory()
    integrator:set_time_disc(ThetaTimeStep(domainDiscT))
    integrator:set_solver(lsolver) -- todo level desc
    return integrator
end

function util.xbraid.create_factory_ConstStepLinearTimeIntegrator()
    integrator = ConstStepLinearTimeIntegratorFactory()
    integrator:set_time_disc(ThetaTimeStep(domainDiscT))
    integrator:set_num_steps(1) -- todo level dependend
    integrator:set_solver(lsolver) -- todo level desc
    return integrator
end

function util.xbraid.create_factory_FixedStepThetaIntegrator()
    integrator = FixedStepThetaIntegratorFactory()
    integrator:set_domain(domainDiscT)
    integrator:set_solver(lsolver) -- todo level desc
    for i = 1, max_level do
        integrator:set_level_num_steps(i-1, orderOrTheta) -- todo level desc
    end
    return integrator
end


function util.xbraid.create_factory_ThetaIntegrator()
    integrator = ThetaIntegratorFactory()
    integrator:set_domain(domainDiscT)
    integrator:set_solver(lsolver) -- todo level desc

    for i = 1, max_level do
        integrator:set_level_theta(i-1, orderOrTheta) -- todo level desc
    end
    return integrator
end


function util.xbraid.create_factory_SimpleIntegrator()
    integrator = SimpleIntegratorFactory()
    integrator:set_domain_disc(domainDiscT)
    integrator:set_solver(nlsolver) -- todo level desc
    integrator:set_dt_min(endTime / 131072) -- todo level desc
    integrator:set_dt_max(endTime)
    integrator:set_reduction_factor(0.2)
    return integrator
end

-- ----------------------------------------------------------------------------
-- Driver
-- ----------------------------------------------------------------------------

function util.xbraid.create_simple_driver_setup()

end


function util.xbraid.create_finecoarse_driver_setup()

end


function util.xbraid.create_finecoarsebase_driver_setup()

end


function util.xbraid.create_leveldependend_driver_setup()

end


-- ----------------------------------------------------------------------------
-- Driver, Integrator, Solver combinations
-- ----------------------------------------------------------------------------
--[[ ResidualStepper -
    m_default_time_step - ThetaTimeStep
    linSolver


ThetaIntegratorNonlinear -
    default_integrator NLThetaIntegrator
    list<integrator> NLThetaIntegrator


IntegratorFactory -
    m_fine_time_integrator_factory - IntegratorFactory
    m_coarse_time_integrator_factory - IntegratorFactory

Integrator -
    default_integrator - ITimeIntegrator
    list<integrator> - ITimeIntegrator

BasicDriver -
    default_integrator - IResidualTimeIntegrator
    list<integrator> IResidualTimeIntegrator
    => Theta Time Step --]]