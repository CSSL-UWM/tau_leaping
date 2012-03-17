#pragma once

#ifdef  NDEBUG

#define my_assert(arg)((void)0)

#else//NDEBUG

#ifndef __CUDA_ARCH__
#include <cassert>
#define my_assert(arg) assert(arg)
#elif __CUDA_ARCH__ >= 20
//use this definitions for nvcc v<4.1
#define my_assert(arg) \
	if(!(arg))\
	{\
		printf("cuda assertion failed at file %s, line %d\n", __FILE__, __LINE__);\
	}
//use this definitions for nvcc v>=4.1
//#define my_assert(arg) assert(arg)
#else//__CUDA_ARCH__<20
#define my_assert(arg) ((void)0)
#endif//__CUDA_ARCH__ ndef

#endif//NDEBUG
