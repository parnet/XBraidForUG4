//todo check
#ifndef UGPLUGIN_XBRAIDFORUG4_CORE_SPACE_TIME_COMMUNICATOR_H
#define UGPLUGIN_XBRAIDFORUG4_CORE_SPACE_TIME_COMMUNICATOR_H

#include <unistd.h>
//2025-03 #include <cassert>
#include <common/assert.h>
#ifdef OPENMP_THREADS
#include <omp.h>
#endif
#include "pcl/pcl_comm_world.h"

//2025-03 #include "libs/braid/braid/braid.hpp"

namespace ug{ namespace xbraid {

    class SpaceTimeCommunicator {
    public:
        //--------------------------------------------------------------------------------------------------------------

        MPI_Comm GLOBAL = PCL_COMM_WORLD;
        MPI_Comm TEMPORAL = PCL_COMM_WORLD;
        MPI_Comm SPATIAL = PCL_COMM_WORLD;

        bool verbose = true;
        int globalsize = 1;
        int temporalsize = 1;
        int spatialsize = 1;


        //--------------------------------------------------------------------------------------------------------------

        SpaceTimeCommunicator() = default;
        virtual ~SpaceTimeCommunicator() = default;

        //--------------------------------------------------------------------------------------------------------------

        void split(int numSpatialProcesses) {

            int world_size, myid;
            MPI_Comm_size(PCL_COMM_WORLD, &world_size);
            UG_ASSERT(world_size % numSpatialProcesses == 0, "process_x * process_t != total_process");

            GLOBAL = PCL_COMM_WORLD;

            globalsize = world_size;
            spatialsize = numSpatialProcesses;
            temporalsize = world_size / numSpatialProcesses;
            //std::cout << "\033[1;32m";
            if (verbose) {
                std::cout << "World size before splitting is:\t" << world_size << std::endl;
            }



            MPI_Comm_rank(GLOBAL, &myid);
            const int xcolor = myid / numSpatialProcesses;
            const int tcolor = myid % numSpatialProcesses;

            MPI_Comm_split(GLOBAL, xcolor, myid, &SPATIAL);
            MPI_Comm_split(GLOBAL, tcolor, myid, &TEMPORAL);

            if (verbose) {
                MPI_Comm_size(GLOBAL, &world_size);
                std::cout << "World size after splitting is:\t" << world_size << std::endl;
                MPI_Comm_size(TEMPORAL, &world_size);
                std::cout << "... with temporal world size:\t" << world_size << std::endl;
                MPI_Comm_size(SPATIAL, &world_size);
                std::cout << "... and spatial world size:\t" << world_size << std::endl << std::endl;
            }
            //std::cout << "\033[0m";

            PCL_COMM_WORLD = SPATIAL; // replaces ugs world communicator with the communicator for spatial
        }

        void set_openmp(int mode, int thread_num) {
#ifdef OPENMP_THREADS
            int this_core_count = 8;
            int this_thread_count = 16;
            if (mode == 0) { // by thread num
                omp_set_num_threads(thread_num);

            } else if (mode == 1 ) { // use full core count
                thread_num =  this_core_count / this->globalsize;

            } else if ( mode == 2 ) { // use full thread cound
                thread_num =  this_thread_count / this->globalsize;
            }
            std::cout << "omp ~~ " << thread_num << std::endl;
#endif

        }

        void set_node_cores() {

        }

        void set_node_threads() {

        }

        void unsplit() {
            PCL_COMM_WORLD = GLOBAL; // reset the world communicator
            SPATIAL = PCL_COMM_WORLD;
            TEMPORAL = PCL_COMM_WORLD;
        }

        int get_global_size() const {
            return globalsize;
        }

        int get_temporal_size() const {
            return temporalsize;
        }

        int get_spatial_size() const {
            return spatialsize;
        }

        int get_temporal_rank() const {
            int rank = 0;
            MPI_Comm_rank(TEMPORAL, &rank);
            return rank;
        }

        int get_spatial_rank() const {
            int rank = 0;
            MPI_Comm_rank(SPATIAL, &rank);
            return rank;
        }

        int get_global_rank() const {
            int rank = 0;
            MPI_Comm_rank(GLOBAL, &rank);
            return rank;
        }

        void sleep(int microseconds) {
            usleep(microseconds);
        }

        //--------------------------------------------------------------------------------------------------------------
    };
}}
#endif