#ifndef UGPLUGIN_XBRAIDFORUG4_INTERFACE_SPATIAL_NORM_H
#define UGPLUGIN_XBRAIDFORUG4_INTERFACE_SPATIAL_NORM_H

#include "lib_disc/function_spaces/grid_function.h"





namespace ug{ namespace xbraid {

/**
 * @brief Interface for computing norms of grid functions.
 *
 * This interface defines the contract for computing various "norms"
 * (e.g. L2, H1, or user-defined) of a GridFunction object.
 *
 * The main purpose of this interface is to provide a uniform, extendable
 * entry point to norm computations, without coupling directly to the
 * core implementation of ug4.
 *
 * Rationale:
 * ----------
 * - ugcore and limex already provides functionality to evaluate norms,
 *   but it is not aware of advanced execution contexts such as
 *   time-parallel methods
 * - By exposing this dedicated interface, plugins can:
 *     * use the same norm computations as the core,
 *     * override or adapt them in contexts where time-dependence matters,
 *     * provide additional diagnostics or logging.
 * - This separation also avoids direct dependencies on the internal
 *   representation of norm evaluation in the core.
 *
 * Notes:
 * ------
 * - This interface could theoretically be replaced by direct use of the
 *   core norm evaluation routines. However, doing so would restrict the
 *   ability of plugins to customize output in parallel-in-time workflows.
 * - Implementers are encouraged to delegate to the core routines whenever
 *   possible to ensure consistency, and only add additional layers when
 *   time-parallel or context-specific behavior is required.
 *
 * Typical usage:
 * --------------
 * - A plugin obtains a GridFunction and passes it to an implementation of
 *   this interface to compute its norm.
 * - For time-parallel plugins, the implementation may incorporate temporal
 *   decomposition, aggregation, or metadata collection in addition to the
 *   core norm computation.
 *
 */

    template <typename TDomain, typename TAlgebra>
    class BraidSpatialNorm {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        //--------------------------------------------------------------------------------------------------------------
    protected:
        BraidSpatialNorm() = default;

    public:
        virtual ~BraidSpatialNorm() = default;

        //--------------------------------------------------------------------------------------------------------------

        virtual double norm(SP_GridFunction u) = 0;

        //--------------------------------------------------------------------------------------------------------------
    };


}}
#endif