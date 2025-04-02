function util.xbraid.create_spatial_grid_transfer(desc,inst)
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
    sgt = SpatialGridTransfer()
    sgt:set_transfer()
    sgt:set_domain(inst.domain_disc)
    sgt:set_approx_space(inst.approx_space)
    if inst.transfer.transfer ~= nil then
        sgt:set_transfer(inst.transfer.transfer)          -- todo userdata, string construction?
    elseif inst.transfer.prolongation ~= nil and inst.transfer.restriction ~= nil then
        sgt:set_prolongation(inst.transfer.prolongation)  -- todo userdata, string construction?
        sgt:set_restriction(inst.transfer.restriction)    -- todo userdata, string construction?
    else
        print("[ ERROR ] transfer or (prolongation and restriction) must be set into *inst*")
        print("util.xbraid.create_spatial_grid_transfer")
        exit()
    end

    return sgt
end
