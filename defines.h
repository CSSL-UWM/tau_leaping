#ifndef _DEFINES_H_
#define _DEFINES_H_

//some usefull defenitions

#include <cusparse.h>

//#define LOG_ENABLED
//#define EXTRA_CHECK
//#define INTERNAL_TIMING

#define DEVICE

#ifdef DEVICE

#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable: 4100 4512 4127)
#endif

#include <thrust/device_vector.h>
	typedef thrust::device_vector<int> vec_int_t;
	typedef thrust::device_vector<char> vec_char_t;
	typedef thrust::device_vector<float> vec_float_t;
	typedef thrust::device_vector<bool> vec_bool_t;
#elif defined HOST
#include <thrust/host_vector.h>
	typedef thrust::host_vector<int> vec_int_t;
	//typedef thrust::host_vector<char> vec_char_t;
	typedef thrust::host_vector<float> vec_float_t;
#endif

#ifdef _MSC_VER
#pragma warning (pop)
#endif


extern cusparseHandle_t g_handle;

namespace std
{
/*#ifdef MSVC
	using namespace tr1;
#endif*/
}

#define CHECK_CUDA_ERRORS() \
do{ \
    cudaError_t cudaErr=cudaGetLastError();\
    if(cudaSuccess!=cudaErr){\
        fprintf(stderr, "CUDA error detected at file: %s, line %d\n" \
                "error code is: %d, error description is: %s\n", \
                __FILE__, __LINE__, cudaErr, cudaGetErrorString(cudaErr));\
        exit(EXIT_FAILURE);\
    }\
    \
}while(0)

#ifndef NDEBUG
#define DEBUG_CHECK_CUDA_ERRORS() CHECK_CUDA_ERRORS()
#else
#define DEBUG_CHECK_CUDA_ERRORS()
#endif


#endif
