#ifndef UGPLUGIN_XBRAIDFORUG4_MATH_VECTOR_HPP
#define UGPLUGIN_XBRAIDFORUG4_MATH_VECTOR_HPP

//2025-03 #include <stddef.h> // size_t
#include <cmath> // log, exp




namespace ug {
    // calculates dest = v2-v1
    template<typename vector_t>
    inline void VectorScaleAdd(vector_t &dest, double alpha, const vector_t &v1, double beta, const vector_t &v2) {
        const size_t dest_size = dest.size();
        for(size_t i=0; i < dest_size; i++) {
            dest[i] = alpha * v1[i] + beta * v2[i];
        }
    }


    template<typename TDomain, typename TAlgebra>
    void VecScaleAdd(GridFunction<TDomain,TAlgebra> &dest, double alpha, const GridFunction<TDomain,TAlgebra> &v1, double beta, const GridFunction<TDomain,TAlgebra> &v2) {
        using TVector = typename GridFunction<TDomain,TAlgebra>::vector_type ;
        VecScaleAdd((TVector&) dest, alpha, (TVector&) v1,beta,(TVector&) v2);
    }
    /*
    // calculates dest = v2-v1
    inline void VecSub(double &dest, const double &v1, const double &v2) {
        dest = v2 - v1;
    }

    // calculates dest = v2-v1
    template<typename vector_t>
    inline void VecSub(vector_t &dest, const vector_t &v1, const vector_t &v2) {
        const size_t dest_size = dest.size();
        for(size_t i=0; i < dest_size; i++)
            VecSub(dest[i], v1[i], v2[i]);
    }
    */
}
#endif