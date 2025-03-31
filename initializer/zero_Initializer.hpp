#ifndef UGPLUGIN_XBRAIDFORUG4_INITIALIZER_ZERO_INITIALIZER_HPP
#define UGPLUGIN_XBRAIDFORUG4_INITIALIZER_ZERO_INITIALIZER_HPP

#include "interface/initializer.hpp"





namespace ug{ namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class ZeroInitializer : public BraidInitializer<TDomain, TAlgebra> {
    public:
        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;

        //--------------------------------------------------------------------------------------------------------------

        ZeroInitializer() = default;

        ~ZeroInitializer() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void initialize(SP_GridFunction& u, number time) override {
            if (time == this->t_start) {
                u = this->m_u0->clone();
            } else {
                u = this->m_u0->clone_without_values();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
    };

}}
#endif