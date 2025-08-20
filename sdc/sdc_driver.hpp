#pragma once


#define __likely__ [[likely]]
#define __unlikely__ [[unlikely]]

namespace ug {
    namespace xbraid {
/*
template <typename TDomain, typename TAlgebra>
class SDC_Object {
    using T_GridFunction = ug::GridFunction<TDomain, TAlgebra>;
    using SP_GridFunction = SmartPtr<T_GridFunction>;

    using T_VectorValueType = typename TAlgebra::vector_type::value_type;

    using SP_SDC_Object = SmartPtr<SDC_Object>;

    SP_GridFunction gf;
    void * send_buffer;
    void * recv_buffer;
    MPI_Request send_request = MPI_REQUEST_NULL;
    MPI_Request recv_request = MPI_REQUEST_NULL;

    SDC_Object();
    SDC_Object(SP_SDC_Object &init);
    ~SDC_Object() {
        if (recv_request != MPI_REQUEST_NULL) {
            int flag = 0;
            MPI_Test(&recv_request, &flag, MPI_STATUS_IGNORE);

            if (!flag) {
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            }

            // Once complete, free old buffer
            free(recv_buffer);
            recv_buffer = nullptr;
            recv_request = MPI_REQUEST_NULL;
        }
    }

    void add(SP_SDC_Object other){
        VecScaleAdd(*gf, 1, *other.gf, 1, *gf);
        }

    void sub(SP_SDC_Object other) {
        VecScaleAdd(*gf, -1, *other.gf, 1, *gf);
    }

    void rmul(double other){
        VecScale(*gf, other);
    }

    void abs() {
        return this->gf->norm(); // or problem related ? set norm?
    }

    int get_buf_size();;

    void pack(void * buffer, T_GridFunction *u_ref, int * buffer_size);

    int buf_pack(void* buffer);

    MPI_Request isend(int dest, int tag, MPI_Comm comm) {
        int buf_size = this->get_buf_size();
        if (send_request != MPI_REQUEST_NULL) {
            int flag = 0;
            MPI_Test(&send_request, &flag, MPI_STATUS_IGNORE);

            if (!flag) {
                MPI_Wait(&send_request, MPI_STATUS_IGNORE);
            }

            // Once complete, free old buffer
            free(send_buffer);
            send_buffer = nullptr;
            send_buffer = MPI_REQUEST_NULL;
        }

        void *send_buf = malloc(buf_size); // todo free space
        buf_pack(send_buf);
        MPI_Request request;
        MPI_Isend(send_buf, buf_size, MPI_BYTE, dest, tag, comm, &request);
        return send_request;
    }

    MPI_Request irecv(int source, int tag, MPI_Comm comm) {
        int buf_size = this->get_buf_size();  // or receive from sender beforehand

        if (recv_request != MPI_REQUEST_NULL) {
            int flag = 0;
            MPI_Test(&recv_request, &flag, MPI_STATUS_IGNORE);

            if (!flag) {
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            }

            // Once complete, free old buffer
            free(recv_buffer);
            recv_buffer = nullptr;
            recv_request = MPI_REQUEST_NULL;
        }

        void *recv_buf = malloc(buf_size);  // todo: store this
        MPI_Irecv(recv_buf, buf_size, MPI_BYTE, source, tag, comm, &recv_request);
        recv_buffer = recv_buf;  // save buffer if needed later

        return recv_request;


    }// communication
    void bcast(SP_SDC_Object other); // communication
};


class SDC_Driver {
public:

    bool initialized_eval = false;
    bool initialized_solve = false;

    Problem();

    void initialize_eval_f(){
    }

    void initialize_solver_system(){
    }

    me_type eval_f(u,t){
        // construct and cache
        if(!initialized_eval) __unlikely__ {
            initialize_eval_f();
        }
        // todo stiffness_matrix * u ( + source_function(t)?)
    }

    me_type solve_system(rhs,dt,u0,t){
        if(!initialized_solve) __unlikely__ {
            initialize_solver_system();
        }
    }
};


template<typename TDomain, typename TAlgebra>
int SDC_Object<TDomain, TAlgebra>::buf_pack(void *buffer) {
#ifdef FEATURE_SPATIAL_REFINE
    __debug(std::cout << "GridFunctionBaseDriver::BufPack" << std::endl);

    /* unpack variables* /
    SP_GridFunction* u_ref = &this->gf;

    int buffer_size = 0;

    auto* chBuffer = static_cast<byte *>(buffer);
    const int spatial_level = u_ref->get()->grid_level().level();
    __debug(std::cout << "Spatial Level: " << spatial_level<< std::endl << std::flush);
    memcpy(chBuffer + buffer_size, &spatial_level, sizeof(int));
    buffer_size += sizeof(int);


    uint mask = u_ref->get()->get_storage_mask();
    __debug(std::cout << "Storage Mask: " << mask<< std::endl << std::flush);
    memcpy(chBuffer + buffer_size, &mask, sizeof(uint));
    buffer_size += sizeof(uint);

    this->pack(buffer, u_ref->get(), &buffer_size);

    __debug(std::cout << "Buffer Size: " << buffer_size << std::endl << std::flush);
    __send_recv_times( std::cout << "Send t=" << timer.get() << std::endl;);
    return 0;

#else
    __debug(std::cout << "GridFunctionBaseDriver::BufPack" << std::endl);
    int buffer_size = 0; // startposition of gridfunction (will be written first) in buffer

    SP_GridFunction* u_ref = &this->gf;

    this->pack(buffer, u_ref->get(), &buffer_size);
    // buffer filled with size of vector and vector

    __send_recv_times( std::cout << "Send t=" << timer.get() << std::endl;);
    return 0;
#endif
}

template<typename TDomain, typename TAlgebra>
int SDC_Object<TDomain, TAlgebra>::get_buf_size() {
    int size_ptr
            __debug(std::cout << "GridFunctionBaseDriver::BufSize" << std::endl);
    size_ptr = 0;
#ifdef FEATURE_SPATIAL_REFINE
    size_ptr =  0

                +sizeof(int)        // spatial-grid-level
                +sizeof(uint)       // parallel storage mask ( undefined, konsistent, unique, additive)
                +sizeof(size_t)     // number of gridfunction-elements
                +sizeof(T_VectorValueType) * (*this->gf).size();  // size of actual vector
#else
        size_ptr =  sizeof(size_t) // number of gridfunction-elements
                     + (sizeof(T_VectorValueType) * (*this->u0).size());
        // size of actual vector

#endif

    __debug(std::cout << "Buffer Size: " << size_ptr << std::endl << std::flush);

    return size_ptr;
}

template<typename TDomain, typename TAlgebra>
void SDC_Object<TDomain, TAlgebra>::pack(void *buffer, T_GridFunction *u_ref, int *buffer_size) {
#ifdef FEATURE_SPATIAL_REFINE
    auto* chBuffer = static_cast<byte *>(buffer);
    const size_t szVector = u_ref->size();
    __debug(std::cout << "num-elem: " << szVector << std::endl);
    memcpy(chBuffer+ *buffer_size, &szVector, sizeof(size_t)); // first value size of vector
    *buffer_size += sizeof(size_t);

    for (size_t i = 0; i < szVector; i++) {
        memcpy(chBuffer + *buffer_size, &(*u_ref)[i], sizeof(T_VectorValueType)); // array sequentially
        *buffer_size += sizeof(T_VectorValueType);
    }

#else

    byte_t* chBuffer = (byte_t*)buffer;

    size_t szVector = u_ref->size();

    memcpy(buffer, &szVector, sizeof(size_t)); // first value size of vector

    *bufferSize += sizeof(size_t);

    for (size_t i = 0; i < szVector; i++) {
        memcpy(chBuffer + *bufferSize, &(*u_ref)[i], sizeof(T_VectorValueType)); // array sequentially
        *bufferSize += sizeof(T_VectorValueType);
    }
#endif
}
*/
    }
}



