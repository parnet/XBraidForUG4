
util = {}
util.xbraid = {}

function util.xbraid.create_start_value_initializer(desc,inst)
    method = GridFunctionInitializer()
    return method
end

function util.xbraid.create_zero_initializer(desc,inst)
    method = ZeroInitializer()
    return method
end

function util.xbraid.create_normal_initializer(desc,inst)
    method = RandomValueInitializer()
    method:set_parameter_normal(desc.initializer.mean, desc.initializer.std)
    return method
end

function util.xbraid.create_uniform_initializer(desc,inst)
    method = RandomValueInitializer()
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