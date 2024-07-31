#ifndef UG_PLUGIN_XBRAIDFORUG4_PARALLEL_IO_GRIDFUNCTION_H
#define UG_PLUGIN_XBRAIDFORUG4_PARALLEL_IO_GRIDFUNCTION_H

#include <fstream>

#include "lib_disc/function_spaces/grid_function.h"

namespace ug { namespace XBraidForUG4 {

    template<typename TDomain, typename TAlgebra>
    class PIOGridFunction {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;
        typedef typename TAlgebra::vector_type::value_type T_VectorValueType;

        //--------------------------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------

        PIOGridFunction() = default;
        ~PIOGridFunction() = default;

        //--------------------------------------------------------------------------------------------------------------
        void write(SP_GridFunction u, const char *path) {
            int mpi_rank;
            MPI_Comm_rank(PCL_COMM_WORLD, &mpi_rank);
            auto *u_ref = u.get();
            std::ofstream outfile;
            size_t szVector = u_ref->size();
            std::stringstream ss;
            if (mpi_rank == 0){
                ss << path << ".gridfunction";
            } else {
                ss << path << "_p" << mpi_rank << ".gridfunction";
            }
            outfile.open(ss.str().c_str(), std::ios::binary | std::ios::out);
            outfile.write((const char *) &szVector, sizeof(size_t));
            // write the value for each gridpoint
            for (size_t i = 0; i < szVector; i++) {
                outfile.write((const char *) &(*u_ref)[i], sizeof(T_VectorValueType));
            }
            outfile.close();
        }

        void read(SP_GridFunction u, const char *path) {
            int mpi_rank;
            MPI_Comm_rank(PCL_COMM_WORLD, &mpi_rank);
            auto *u_ref = u.get();
            std::ifstream infile;
            std::stringstream ss;
            if (mpi_rank == 0){
                ss << path << ".gridfunction";
            } else {
                ss << path << "_p" << mpi_rank << ".gridfunction";
            }
            infile.open(ss.str().c_str(), std::ios::binary | std::ios::in);
            if (!infile.good()){
                std::cout << "filename: " << path << std::endl << std::flush;
                UG_THROW("file does not exist");
            }
            size_t szVector = 0;
            infile.read((char *) &szVector, sizeof(size_t));
            for (size_t i = 0; i < szVector; i++) {
                T_VectorValueType val = 0;
                infile.read((char *) &val, sizeof(T_VectorValueType));
                (*u_ref)[i] = val;
            }
            infile.close();
        }
        //--------------------------------------------------------------------------------------------------------------
    };

}}

#endif //UG_PLUGIN_XBRAIDFORUG4_PARALLEL_IO_GRIDFUNCTION_H
