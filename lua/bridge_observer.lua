
util = {}
util.xbraid = {}

function util.xbraid.create_time_integrator_observer_collector(desc,inst)
    method = TimeIntegratorObserverCollector()
    method:attach_observer()
    return method
end

function util.xbraid.create_xBraid_time_integrator_observer_collector(desc,inst)
    method = XBraidTimeIntegratorObserverCollector()
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
    method = VTK_Observer(inst.vtk, desc.filename)
    return method
end

function util.xbraid.create_matlab_observer(desc,inst)
    method = MATLAB_Observer()
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
    method = VTK_ProcessObserver(inst.vtk, desc.filename)
    return method
end


function util.xbraid.create_eval_observer(desc,inst)
    method = EvalObserver()
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
    if desc.name == "EvalObserver" then
        observer = util.xbraid.create_eval_observer(desc,inst)
    elseif desc.name == "VTKProcessObserver" then
        observer = util.xbraid.create_vtk_process_observer(desc,inst)
    elseif desc.name == "MATLABObserver" then
        observer = util.xbraid.create_matlab_observer(desc,inst)
    elseif desc.name == "VTKObserver" then
        observer = util.xbraid.create_vtk_observer(desc,inst)
    end
end

-- todo split process observer and observer
function util.xbraid.create_observer(desc,inst)
    if type(desc.observer) == "userdata" then
        inst.observer = desc.observer
    elseif type(desc.observer) == "table" then
        if desc.observer.name == "ProcessObserverCollector" then
            observer = util.xbraid.create_xBraid_time_integrator_observer_collector(desc,inst)
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
            observer = util.xbraid.create_time_integrator_observer_collector(desc,inst)
            inst.observer = observer
            for k, v in pairs(desc.observer.observers) do
                observer = util.xbraid.create_single_observer(v)
                observer.attach_observer(observer)
            end
        else
            observer = util.xbraid.create_single_observer(desc.observer.name)
            inst.observer = observer
        end
    end
end