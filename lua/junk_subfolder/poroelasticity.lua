util.poro = util.poro or {}



function util.poro.create_conv_check(approxSpace, dim)
    -- todo ø
    local p0 = 1.0 -- scaling factor for every component

    local cmpConvCheck = CompositeConvCheck(approxSpace)
    cmpConvCheck:set_component_check("ux", p0 * 1e-14, 1e-6)
    cmpConvCheck:set_component_check("uy", p0 * 1e-14, 1e-6)
    if (dim == 3) then
        cmpConvCheck:set_component_check("uz", p0 * 1e-14, 1e-6)
    end
    cmpConvCheck:set_component_check("p", p0 * 1e-14, 1e-6)


    cmpConvCheck:set_maximum_steps(100)
    cmpConvCheck:set_verbose(true)
    cmpConvCheck:set_supress_unsuccessful(true)

    return cmpConvCheck
end



function util.poro.create_smoother(smootherID)
    -- todo ø
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



function util.poro.create_uzawa(sSchurCmp, aiForward, aiSchur, aiBackward, uzawaSchurUpdateOp, uzawaSchurWeight)
    -- todo ø
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