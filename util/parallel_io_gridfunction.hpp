#ifndef UGPLUGIN_XBRAIDFORUG4_UTIL_PARALLEL_IO_GRIDFUNCTION_HPP
#define UGPLUGIN_XBRAIDFORUG4_UTIL_PARALLEL_IO_GRIDFUNCTION_HPP

#include <fstream>
#include "lib_disc/function_spaces/grid_function.h"




namespace ug { namespace xbraid {
    /* note that this class may not be system independent */
    template<typename TDomain, typename TAlgebra>
    class PIOGridFunction {
    public:

        //--------------------------------------------------------------------------------------------------------------

        using T_GridFunction = GridFunction<TDomain, TAlgebra> ;
        using SP_GridFunction = SmartPtr<T_GridFunction> ;
        using T_VectorValueType =  typename TAlgebra::vector_type::value_type ;

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
            outfile.write(reinterpret_cast<const char *>(&szVector), sizeof(size_t));

            // write the value of each gridpoint
            for (size_t i = 0; i < szVector; i++) {
                outfile.write(reinterpret_cast<const char *>(&(*u_ref)[i]), sizeof(T_VectorValueType));
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
                std::cout << "filename: " << ss.str().c_str() << std::endl << std::flush;
                UG_THROW("file does not exist");
            }
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