#ifndef UGPLUGIN_XBRAIDFORUG4_UTIL_IO_GRIDFUNCTION_HPP
#define UGPLUGIN_XBRAIDFORUG4_UTIL_IO_GRIDFUNCTION_HPP

#include <fstream>

#include "lib_disc/function_spaces/grid_function.h"



namespace ug { namespace xbraid {

    template <typename TDomain, typename TAlgebra>
    class IOGridFunction {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SPGridFunction = SmartPtr<T_GridFunction> ;
        using T_VectorValueType = typename TAlgebra::vector_type::value_type ;

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
            outfile.write(reinterpret_cast<const char *>(&szVector), sizeof(size_t));
            for (size_t i = 0; i < szVector; i++) {
                outfile.write(reinterpret_cast<const char *>(&(*u_ref)[i]), sizeof(T_VectorValueType));
            }
            outfile.close();
        }

        void read(SPGridFunction u, const char* path) {
            auto* u_ref = u.get();
            std::ifstream infile;
            infile.open(path, std::ios::binary | std::ios::in);
            //if (!infile.good()) {
            //    std::cout << "filename: " << path << std::endl << std::flush;
            //    UG_THROW("file does not exist");
            //}
            size_t szVector = 0;
            infile.read(reinterpret_cast<char *>(&szVector), sizeof(size_t));
            for (size_t i = 0; i < szVector; i++) {
                T_VectorValueType val = T_VectorValueType(0);
                infile.read(reinterpret_cast<char *>(&val), sizeof(T_VectorValueType));
                (*u_ref)[i] = val;
            }
            infile.close();
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif