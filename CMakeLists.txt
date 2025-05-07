cmake_minimum_required(VERSION 3.5)

set(pluginName XBraidForUG4)
project(UG_PLUGIN_${pluginName})


# ----------  ----------  ---------- ---------- Compiling libraries ---------- ----------  ----------  ----------

# [[ create braid as extra library]]  >>>>>>>>>>>>>>
#if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/libbraid.a")
#    message("Compiling XBRAID")
#    add_custom_command(OUTPUT braidlib
#            COMMAND make braid
#            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid
#
#            COMMAND ${CMAKE_COMMAND} -E copy
#                "${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a"
#                "${CMAKE_CURRENT_BINARY_DIR}/libbraid.a"
#
#            COMMAND make clean
#            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid
#    )
#else ()
#    message(STATUS "XBraidForUG4 - libbraid already exists")
#    add_custom_command(
#            OUTPUT braidlib
#            COMMAND echo :
#    )
#endif ()
# add_custom_target(run ALL DEPENDS braidlib)
# [[ create braid as extra library]]  <<<<<<<<<<<<<<<<<<<<<<

#if(${LIMEX})
#elseif ()
#    message(ERROR "Activate Limex to use ${pluginName}")
#endif ()



if(NOT UG_ROOT_CMAKE_PATH)
    set(UG_ROOT_CMAKE_PATH /storage/pro/ug/ugcore/cmake)
endif(NOT UG_ROOT_CMAKE_PATH)


#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)



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
)


set(SOURCES
        bridge/bridge.cpp
        bridge/transfer.cpp

        tests/class_reusage.cpp

        util/memory_observer.cpp
        util/world_memory.cpp
)

set(HEADERS
        config/pragma.hpp
        config/compile_settings.hpp
        config/execution_settings.hpp

        convergence/adaptive_conv_check.hpp
        convergence/fine_coarse_conv_check.hpp
        convergence/simple_conv_check.hpp
        convergence/spatial_accuracy.hpp

        core/braid_settings.hpp
        core/braid_executor.hpp
        core/braid_vector_struct.hpp
        core/space_time_communicator.hpp
        core/replace_standard_stream.hpp
        core/spatial_accuracy.hpp
        core/discard_standard_stream.hpp

        driver/basic_driver.hpp
        driver/braid_integrator.hpp
        driver/braid_integrator_factory.hpp
        driver/braid_residual_stepper.hpp
        driver/braid_theta_integrator_nonlinear.hpp
        driver/gridfunction_base.hpp
        driver/script_writer_driver.hpp

        factory/bdf_integrator_factory.hpp
        factory/const_step_linear_time_integrator_factory.hpp
        factory/fixed_step_theta_integrator_factory.hpp
        factory/limex_integrator_factory.hpp
        factory/linear_time_integrator_factory.hpp
        factory/simple_integrator_factory.hpp
        factory/theta_integrator_factory.hpp

        initializer/lua_initializer.hpp
        initializer/random_value_initializer.hpp
        initializer/start_value_initializer.hpp
        initializer/zero_Initializer.hpp

        integrator/bdf_integrator.hpp
        integrator/bdf_integrator_nl.hpp
        integrator/experimental_timestep.hpp
        integrator/theta_conststep_integrator.hpp
        integrator/theta_conststep_integrator_nl.hpp
        integrator/theta_integrator.hpp
        integrator/theta_integrator_nl.hpp
        integrator/theta_single_timestep.hpp

        interface/initializer.hpp
        interface/integrator_factory.hpp
        interface/observer_xbraid.hpp
        interface/prolongate.hpp
        interface/residual_timeintegrator.hpp
        interface/restrict.hpp
        interface/spatial_norm.hpp

        math/vector.hpp

        observer/collector_observer.hpp
        observer/io_observer.hpp
        observer/io_process_observer.hpp
        observer/matlab_observer.hpp
        observer/no_observer.hpp
        observer/userdata_process_evaluate_observer.hpp
        observer/vtk_observer.hpp
        observer/vtk_process_observer.hpp
        observer/xb_collector_observer.hpp

        problemset/problemsets.hpp

        spatial_norm/euclidian_norm.hpp

        tests/class_reusage.hpp

        transfer/__package.hpp
        transfer/prolongate.hpp
        transfer/prolongate_elementwise.hpp
        transfer/prolongate_p1.hpp
        transfer/real_coarsen.hpp
        transfer/restrict.hpp
        transfer/restrict_elementwise.hpp
        transfer/restrict_p1.hpp
        transfer/spatial_grid_transfer.hpp

        util/braid_timer.hpp
        util/braid_usage_timer.hpp
        util/debug.hpp
        util/io_gridfunction.hpp
        util/memory_observer.hpp
        util/parallel_io_gridfunction.hpp
        util/parallel_logger.hpp
        util/pragma.hpp
        util/world_memory.hpp
    )
# ----------------------------------------- ----------------------------------------------------------------------------
if(DEFINED ConvectionDiffusion)
    if (${Poroelasticity})
        message(STATUS "[X] XBraid Convection Diffusion")

        set(SOURCES_CONVECTION_DIFFUSION
        )

        set(HEADERS_CONVECTION_DIFFUSION

        )

        set(SOURCES ${SOURCES} ${SOURCES_CONVECTION_DIFFUSION})

    elseif ()
        message(STATUS "[X] XBraid Convection Diffusion")
    endif ()
elseif ()
    message(STATUS "[X] XBraid Convection Diffusion")
endif ()


if(DEFINED Poroelasticity)
    if (${Poroelasticity})
        set(XBraidPoroelasticity ON)
        add_compile_definitions(XBraidPoroelasticity=1)
        message(STATUS "[X] XBraid Poroelasticity")

        set(SOURCES_POROELASTICITY
        )

        set(HEADERS_POROELASTICITY
                problemset/poroelasticity/bridge.h
                problemset/poroelasticity/braid_biot_estimator.h
                problemset/poroelasticity/braid_biot_control.h
                problemset/poroelasticity/braid_biot_precomputed.h
                problemset/poroelasticity/braid_heat_check.h
                problemset/poroelasticity/biot_error_data.h
                problemset/poroelasticity/biot_braid_displacement_norm.h
        )

        set(SOURCES ${SOURCES} ${SOURCES_POROELASTICITY}
                lua/header.hpp)

    else()
        message(STATUS "[ ] XBraidPoroelasticity")
    endif ()
else ()
    message(STATUS "[ ] XBraidPoroelasticity")
endif ()


# ---------------------------------------------------------------------------------------------------------------------


if(NOT UG_ROOT_CMAKE_PATH)
    set(UG_ROOT_CMAKE_PATH /storage/pro/ug/ugcore/cmake)
endif(NOT UG_ROOT_CMAKE_PATH)

include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


if (buildEmbeddedPlugins)
    exportsources(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
    exportsources(${CMAKE_CURRENT_SOURCE_DIR} ${BRAIDSOURCES})
else (buildEmbeddedPlugins)
    add_library(${pluginName} SHARED ${SOURCES} ${BRAIDSOURCES})
    include_directories(${pluginName} PUBLIC "${CMAKE_SOURCE_DIR}/plugins/${pluginName}")

    #add_library(${pluginName} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a)
    # target_link_libraries(${pluginName} ${CMAKE_CURRENT_SOURCE_DIR}/libs/braid/braid/libbraid.a)
    #target_link_libraries(${pluginName} "${CMAKE_CURRENT_BINARY_DIR}/libbraid.a")

    target_link_libraries(${pluginName} ug4)
    target_link_libraries(${pluginName} Limex)
    #add_dependencies(${pluginName} UG_PLUGIN_ConvectionDiffusion)
    #LINK_DIRECTORIES(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    #MESSAGE(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    # target_link_libraries(${pluginName} ConvectionDiffusion)
    if(XBraidPoroelasticity)
        target_link_libraries(${pluginName} Poroelasticity)
    endif ()

    target_include_directories(${pluginName}  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..) # Make plugin folder available
endif (buildEmbeddedPlugins)



