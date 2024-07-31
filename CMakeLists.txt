cmake_minimum_required(VERSION 3.5)

set(pluginName XBraidForUG4)
project(UG_PLUGIN_${pluginName})

# ----------  ----------  ---------- ---------- Compiling libraries ---------- ----------  ----------  ----------
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/libbraid.a")
    MESSAGE("Compiling XBRAID")
    add_custom_command(OUTPUT braidlib
            COMMAND make braid
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid

            COMMAND ${CMAKE_COMMAND} -E copy
                "${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a"
                "${CMAKE_CURRENT_BINARY_DIR}/libbraid.a"

            COMMAND make clean
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid
    )
else ()
    MESSAGE("XBraidForUG4 - libbraid already exists")
    add_custom_command(
            OUTPUT braidlib
            COMMAND echo :
    )
endif ()


add_custom_target(run ALL DEPENDS braidlib)

if(NOT UG_ROOT_CMAKE_PATH)
    set(UG_ROOT_CMAKE_PATH /storage/pro/ug/ugcore/cmake)
endif(NOT UG_ROOT_CMAKE_PATH)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_INIT} -fPIC")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


set(SOURCES

        src/xbraidforug4_plugin.cpp



        src/core/braid_settings.h
        src/core/braid_executor.h
        src/core/braid_vector_struct.h
        src/core/space_time_communicator.h

        src/driver/braid_integrator.h
        src/driver/gridfunction_base.h
        src/driver/basic_driver.h
        src/driver/braid_integrator_factory.h
        src/driver/braid_residual_stepper.h
        src/driver/braid_theta_integrator_nonlinear.h

        src/factory/bdf_integrator_factory.h
        src/factory/const_step_linear_time_integrator_factory.h
        src/factory/fixed_step_theta_integrator_factory.h
        src/factory/limex_integrator_factory.h
        src/factory/linear_time_integrator_factory.h
        src/factory/simple_integrator_factory.h
        src/factory/theta_integrator_factory.h

        src/initializer/start_value_initializer.h
        src/initializer/zero_Initializer.h

        src/integrator/theta_single_timestep.h # B. Sc. Thesis TimeStepper for Heat Equation
        src/integrator/theta_integrator.h
        src/integrator/bdf_integrator.h
        src/integrator/theta_conststep_integrator_nl.h
        src/integrator/bdf_integrator_nl.h
        src/integrator/theta_integrator_nl.h
        src/integrator/theta_conststep_integrator.h



        src/interface/integrator_factory.h
        src/interface/initializer.h
        src/interface/residual_timeintegrator.h
        src/interface/spatial_norm.h

        src/spatial_norm/euclidian_norm.h

        src/util/braid_timer.h
        src/util/braid_usage_timer.h
        src/util/io_gridfunction.h
        src/util/memory_observer.h
        src/util/memory_observer.cpp
        src/util/parallel_io_gridfunction.h
        src/util/parallel_logger.h
        src/util/world_memory.h
        src/util/world_memory.cpp


        src/interface/observer_xbraid.h
        src/observer/userdata_process_evaluate_observer.h
        src/observer/matlab_observer.h
        src/observer/collector_observer.h
        src/common/pragma.h
        src/core/replace_standard_stream.h
    )

set(BRAIDSOURCES
        libs/braid/braid/_braid.h
        libs/braid/braid/access.c
        libs/braid/braid/adjoint.c
        libs/braid/braid/adjoint.h
        libs/braid/braid/base.c
        libs/braid/braid/base.h
        libs/braid/braid/braid.c
        libs/braid/braid/braid.h
        libs/braid/braid/braid.hpp
        libs/braid/braid/braid_defs.h
        # libs/braid/braid/braid_F90_iface.c
        libs/braid/braid/braid_status.c
        libs/braid/braid/braid_status.h
        libs/braid/braid/braid_test.c
        libs/braid/braid/braid_test.h
        libs/braid/braid/communication.c
        libs/braid/braid/distribution.c
        libs/braid/braid/drive.c
        libs/braid/braid/grid.c
        libs/braid/braid/hierarchy.c
        libs/braid/braid/interp.c
        libs/braid/braid/mpistubs.c
        libs/braid/braid/mpistubs.h
        libs/braid/braid/norm.c
        libs/braid/braid/refine.c
        libs/braid/braid/relax.c
        libs/braid/braid/residual.c
        libs/braid/braid/restrict.c
        libs/braid/braid/space.c
        libs/braid/braid/status.h
        libs/braid/braid/step.c
        libs/braid/braid/tape.c
        libs/braid/braid/tape.h
        libs/braid/braid/delta.c
        libs/braid/braid/delta.h
        libs/braid/braid/util.c
        libs/braid/braid/util.h
        libs/braid/braid/uvector.c
        src/core/replace_standard_stream.h
        src/observer/mock_observer.h
        src/observer/vtk_observer.h
        src/observer/io_observer.h
        src/observer/vtk_process_observer.h
        src/observer/io_process_observer.h
        src/observer/xb_collector_observer.h
        src/common/XObject.h
        src/common/XObject.h
)


if(NOT UG_ROOT_CMAKE_PATH)
    set(UG_ROOT_CMAKE_PATH /storage/pro/ug/ugcore/cmake)
endif(NOT UG_ROOT_CMAKE_PATH)

include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if (buildEmbeddedPlugins)
    EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
    EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${BRAIDSOURCES})
else (buildEmbeddedPlugins)
    add_library(${pluginName} SHARED ${SOURCES} ${BRAIDSOURCES})
    #add_library(${pluginName} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a)
    # target_link_libraries(${pluginName} ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a)
    target_link_libraries(${pluginName} "${CMAKE_CURRENT_BINARY_DIR}/libbraid.a")
    target_link_libraries(${pluginName} ug4)
    target_link_libraries(${pluginName} Limex)
    #add_dependencies(${pluginName} UG_PLUGIN_ConvectionDiffusion)
    LINK_DIRECTORIES(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    #MESSAGE(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    # target_link_libraries(${pluginName} ConvectionDiffusion)
    
endif (buildEmbeddedPlugins)



