#include "world_memory.hpp"

#include <vector>
#include <mpi.h>
#include "core/space_time_communicator.hpp"
//2025-03 #include "lib_disc/function_spaces/grid_function.h"

#include "memory_observer.hpp"


namespace ug { namespace xbraid {

        unsigned long get_world_memory_consumed() {
            unsigned long world_memory = 0;
            const auto this_proc_memory = get_physical_memory_consumed();
            MPI_Allreduce(&this_proc_memory,
                          &world_memory,
                          1, MPI_UNSIGNED_LONG, MPI_SUM,
                          MPI_COMM_WORLD);
            return world_memory;
        }

        unsigned long get_world_memory_peak(){
            unsigned long world_memory = 0;
            const auto this_proc_memory = get_physical_memory_peak();

            MPI_Allreduce(&this_proc_memory,
                          &world_memory,
                          1, MPI_UNSIGNED_LONG, MPI_SUM,
                          MPI_COMM_WORLD);

            return world_memory;
        }

        unsigned long get_spatial_memory_consumed(){
            unsigned long world_memory = 0;
            const auto this_proc_memory = get_physical_memory_consumed();

            MPI_Allreduce(&this_proc_memory,
                          &world_memory,
                          1, MPI_UNSIGNED_LONG, MPI_SUM,
                          PCL_COMM_WORLD);

            return world_memory;
        }

        void get_world_memory_distribution(){
            int size = 0;
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            auto distribution = std::vector<unsigned long>();
            distribution.resize(size);
            auto* sbuf = static_cast<unsigned long *>(malloc(1 * sizeof(unsigned long)));
            auto* rbuf = static_cast<unsigned long *>(malloc(size * sizeof(unsigned long)));
            sbuf[0] = get_physical_memory_consumed();
            MPI_Allgather(sbuf, 1, MPI_UNSIGNED_LONG,
                          rbuf, 1, MPI_UNSIGNED_LONG,
                          MPI_COMM_WORLD);
            unsigned long summ = 0.0;
            for (int i = 0; i < size; i++) {
                summ += rbuf[i];
                std::cout << i << ": " << rbuf[i] << std::endl;
            }
            std::cout << "Total : " << summ << std::endl;
        }

        void get_spatial_memory_distribution()   {
            int size = 0;
            MPI_Comm_size(PCL_COMM_WORLD, &size);
            auto distribution = std::vector<unsigned long>();
            distribution.resize(size);
            auto* sbuf = static_cast<unsigned long *>(malloc(1 * sizeof(unsigned long)));
            auto* rbuf = static_cast<unsigned long *>(malloc(size * sizeof(unsigned long)));
            sbuf[0] = get_physical_memory_consumed();
            MPI_Allgather(sbuf, 1, MPI_UNSIGNED_LONG,
                          rbuf, 1, MPI_UNSIGNED_LONG,
                          PCL_COMM_WORLD);
            unsigned long summ = 0.0;
            for (int i = 0; i < size; i++) {
                summ += rbuf[i];
                std::cout << i << ": " << rbuf[i] << std::endl;
            }
            std::cout << "Total : " << summ << std::endl;
        }

    }}
