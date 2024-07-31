//todo check
#ifndef UG_PLUGIN_XBRAIDFORUG4_SPACE_TIME_COMMUNICATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_SPACE_TIME_COMMUNICATOR_H

#include <unistd.h>
#include <cassert>
#include <common/assert.h>

#include "pcl/pcl_comm_world.h"

#include "../../libs/braid/braid/braid.hpp"



namespace ug { namespace XBraidForUG4 {

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

            PCL_COMM_WORLD = SPATIAL; // replaces ugs world communicator with the communicator for spatial
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
#endif //UG_PLUGIN_XBRAIDFORUG4_SPACE_TIME_COMMUNICATOR_H
