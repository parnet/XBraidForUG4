#include "world_memory.h"

#include <vector>

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"

#include "memory_observer.h"



namespace ug { namespace XBraidForUG4 {

        unsigned long get_world_memory_consumed() {
            unsigned long world_memory = 0;
            auto this_proc_memory = get_physical_memory_consumed();
            MPI_Allreduce(&this_proc_memory,
                          &world_memory,
                          1, MPI_UNSIGNED_LONG, MPI_SUM,
                          MPI_COMM_WORLD);
            return world_memory;
        }

        unsigned long get_world_memory_peak(){
            unsigned long world_memory = 0;
            auto this_proc_memory = get_physical_memory_peak();

            MPI_Allreduce(&this_proc_memory,
                          &world_memory,
                          1, MPI_UNSIGNED_LONG, MPI_SUM,
                          MPI_COMM_WORLD);

            return world_memory;
        }

        unsigned long get_spatial_memory_consumed(){
            unsigned long world_memory = 0;
            auto this_proc_memory = get_physical_memory_consumed();

            MPI_Allreduce(&this_proc_memory,
                          &world_memory,
                          1, MPI_UNSIGNED_LONG, MPI_SUM,
                          PCL_COMM_WORLD);

            return world_memory;
        }

        void get_world_memory_distribution(){
            int size = 0;
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            std::vector<unsigned long> distribution = std::vector<unsigned long>();
            distribution.resize(size);
            unsigned long* sbuf = (unsigned long*)malloc(1 * sizeof(unsigned long));
            unsigned long* rbuf = (unsigned long*)malloc(size * sizeof(unsigned long));
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
            std::vector<unsigned long> distribution = std::vector<unsigned long>();
            distribution.resize(size);
            unsigned long* sbuf = (unsigned long*)malloc(1 * sizeof(unsigned long));
            unsigned long* rbuf = (unsigned long*)malloc(size * sizeof(unsigned long));
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
