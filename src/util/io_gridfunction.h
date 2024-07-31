#ifndef UG_PLUGIN_XBRAIDFORUG4_IOGRIDFUNCTION_H
#define UG_PLUGIN_XBRAIDFORUG4_IOGRIDFUNCTION_H

#include <fstream>

#include "lib_disc/function_spaces/grid_function.h"


namespace ug { namespace XBraidForUG4 {

    template <typename TDomain, typename TAlgebra>
    class IOGridFunction {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SPGridFunction;
        typedef typename TAlgebra::vector_type::value_type T_VectorValueType;

        //--------------------------------------------------------------------------------------------------------------

        IOGridFunction() = default;
        ~IOGridFunction() = default;

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        void write(SPGridFunction u, const char* path) {
            auto* u_ref = u.get();
            std::ofstream outfile;
            size_t szVector = u_ref->size();

            outfile.open(path, std::ios::binary | std::ios::out);
            outfile.write((const char*)&szVector, sizeof(size_t));
            for (size_t i = 0; i < szVector; i++) {
                outfile.write((const char*)&(*u_ref)[i], sizeof(T_VectorValueType));
            }
            outfile.close();
        }

        void read(SPGridFunction u, const char* path) {
            auto* u_ref = u.get();
            std::ifstream infile;
            infile.open(path, std::ios::binary | std::ios::in);
            if (!infile.good()) {
                std::cout << "filename: " << path << std::endl << std::flush;
                UG_THROW("file does not exist");
            }
            size_t szVector = 0;
            infile.read((char*)&szVector, sizeof(size_t));
            for (size_t i = 0; i < szVector; i++) {
                T_VectorValueType val = 0;
                infile.read((char*)&val, sizeof(T_VectorValueType));
                (*u_ref)[i] = val;
            }
            infile.close();
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif //UG_PLUGIN_XBRAIDFORUG4_IOGRIDFUNCTION_H
