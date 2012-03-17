#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

//this file contains a set of common functions

#include <vector>
#include <stdio.h>
#include "realtype.h"

#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable: 4100 4512 4127)
#endif

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#ifdef _MSC_VER
#pragma warning (pop)
#endif


#include <cmath>
#include "d_rng.cuh"
#include "my_assert.h"

__host__ __device__
inline int DivUp(int i1, int i2)
{
	my_assert(i1>0);
	return(i1-1)/i2+1;
}


//Box-Muller algorithm for creating a normal-distributed value
//out of two uniform-distributed
inline __host__ __device__
float boxMuller ( float mean, float variance, uint4 &seed1, uint4 &seed2)
{
	float s, u, v;
    do {
	   float a=HybridTausRng(&seed1);
	   float b=HybridTausRng(&seed2);
                   // a,b uniform in (0,1)
       u = 2*a-1;
       v = 2*b-1;  // u,v uniform in (-1,1)
       s = u*u + v*v;
    } while (s >= 1);
  
                   // (u,v) uniform in unit circle 
    float  x = u*sqrt(-2*log(s)/s); 
   // float y = v*sqrt(-2*log(s)/s); 
    return ( mean + sqrt(variance)*x );
}

//get random values in seeds array
bool FillRandoms(char* seeds, unsigned int size);

//read a file into memory (for XML parser)
void FileToMem(const char *Name, std::vector<char> &buffer);

//check if a file exists
bool FileExists(const char *fname);

//printing out a device vector for debugging
template <typename T, typename OutT>
void print_vector(const thrust::device_vector<T> &v, const char *delim=" ")
{
	thrust::host_vector<OutT> h=v;
	for(size_t i=0; i<h.size(); ++i)
	{
		std::cout<<h[i];
		if(i!=h.size()-1)
			std::cout<<delim;
	}
}

//compute propensity out of reaction's constant, reaction's type, 
inline __device__ __host__
PropReal_t Propensity(float k, char type, int x0, int x1)
{
	//return k*x0*(x1*(type<<2)+(type<<1)+0.5f*(2*x0-1)*type);
	switch (type)
	{
	case 3:
		return k*x0*(x0-1)/2.f;
	case 2:
		return k*x0*x1;
	case 1:
		return k*x0;
	}
#ifndef NDEBUG
	printf("Wrong type occured: %d\n", type);
	my_assert(false);
#endif
	return -1.f;
}

//compute if a the reaction is critical or not
inline __device__ __host__
void IsCriticalAndPropensity(char type, int x0, int x1, int nc, float k, 
							 float *prop, bool *isCrit)
{
	switch(type)
	{
	case 1:
		*isCrit=x0<nc;
		*prop=k*x0;
		return;
	case 2:
		*isCrit=x0<nc||x1<nc;
		*prop=k*x0*x1;
		return;
	case 3:
		*isCrit=x0*2<nc;
		*prop=k*x0*(x0-1)/2.f;
		return;
	default:
		my_assert(false);
		return;
	}

}

//pack reactant's constant and column index into one integer
inline __host__ __device__
int pack(char v, int i)
{
	return (i<<3)+((char)(v+2)&7);
}

template <typename SmallPart_t, int bits>
__forceinline__ __host__ __device__
void split_templ(int vci, SmallPart_t &v, int &ci)
{
	ci=vci>>bits;
	v=SmallPart_t(vci&((1<<bits)-1));
}

//split an integer to reactant's constant and column index
__forceinline__ __host__ __device__
void split(int vci, char &v, int &ci)
{
	split_templ<char, 3>(vci, v, ci);
	v-=2;
}

class ValueAndInd2Int
{
public: 
	__host__
		ValueAndInd2Int(){}
	inline __host__ __device__
	int operator()(char fa, int sa)const
	{
		
		return pack(fa, sa);
	}
};

__device__ __host__
inline
int log2(int val)
{
	my_assert(val!=0);
	int targetlevel = 0;
	while (val >>= 1) ++targetlevel;
	return targetlevel;
};


__device__ __host__
inline
int nextPowerOf2(int val)
{
	int lg2=log2(val);
	int p2=1<<lg2;
	if(p2<val)
		p2*=2;

	return p2;
}

#endif