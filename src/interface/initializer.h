#ifndef UG_PLUGIN_XBRAIDFORUG4_INITIALIZER_H
#define UG_PLUGIN_XBRAIDFORUG4_INITIALIZER_H

namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class BraidInitializer {
    public:
        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        SP_GridFunction m_u0 = SPNULL; // u0 and t_start for star value problem
        double t_start = 0.0;

        //--------------------------------------------------------------------------------------------------------------

        BraidInitializer() = default;

        virtual ~BraidInitializer() = default;

        //--------------------------------------------------------------------------------------------------------------
        virtual void initialize(SP_GridFunction& u, number time) = 0;

        virtual void set_start_values(SP_GridFunction u, number time) {
            this->m_u0 = u;
            this->t_start = time;
        }

        //--------------------------------------------------------------------------------------------------------------
    };

}}

#endif //UG_PLUGIN_XBRAIDFORUG4_INITIALIZER_H
