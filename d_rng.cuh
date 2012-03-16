#ifndef __D_RNG_CUH__
#define __D_RNG_CUH__

//hybrid taus random number generator

#define NEXT_RAND(x)         ((x) * 1664525 + 1013904223UL)

typedef unsigned int uint;

__host__        uint    CreateRandomSeedsUint4(uint4** seeds, uint size);

__device__ __host__ __forceinline__ uint  LCGStep(unsigned int* seed)
{
       *seed = NEXT_RAND(*seed);
       return *seed;
}

__device__ __host__ __forceinline__ uint  TausStep(uint* z, int S1, int S2, int S3, uint M)
{
       uint    b = ((((*z) << S1) ^ (*z)) >> S2);

       *z = ((((*z) & M) << S3) ^ b);

       return  *z;
}

__device__ __host__ inline float   HybridTausRng(uint4*    state)
{
       return 2.3283064365387e-10f * (
                       TausStep(&state->x, 13, 19, 12, 4294967294UL) ^
                       TausStep(&state->y, 2, 25, 4, 4294967288UL) ^
                       TausStep(&state->z, 3, 11, 17, 4294967280UL) ^
                       LCGStep(&state->w)
       );
}

__device__ __host__ __forceinline__ float   HybridTausRngReduced(uint4*    state)
{
       return 2.3283063e-10f * (
                       TausStep(&state->x, 13, 19, 12, 4294967294UL) ^
                       TausStep(&state->y, 2, 25, 4, 4294967288UL) ^
                       TausStep(&state->z, 3, 11, 17, 4294967280UL) ^
                       LCGStep(&state->w)
       );
}

__device__ inline uint HybridTausRngInt(uint4* state)
{
       return (
                   TausStep(&state->x, 13, 19, 12, 4294967294UL) ^
                   TausStep(&state->y, 2, 25, 4, 4294967288UL) ^
                       TausStep(&state->z, 3, 11, 17, 4294967280UL) ^
                       LCGStep(&state->w)
       );
}

#endif
