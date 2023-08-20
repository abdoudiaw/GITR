#ifndef _SURFACE_
#define _SURFACE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include "Surfaces.h"
#include <cmath>
#include "pusher.h"

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#endif

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

CUDA_CALLABLE_MEMBER
gitr_precision interp2d(gitr_precision x, gitr_precision z, int nx, int nz,
                        std::vector<double>& gridx, const std::vector<double>& gridz, std::vector<double>& data);
gitr_precision interp2d(gitr_precision x, gitr_precision z, int nx, int nz,
                        const std::vector<double>& gridx, const std::vector<double>& gridz, const std::vector<double>& data);
struct surfaceReactions {
    Particles * particles;
    int nLines;
    Boundary * boundaryVector;
    Surfaces * surfaces;
#if __CUDACC__
        curandState *state;
#else
        std::mt19937 *state;
#endif
    Flags * flags;

    surfaceReactions(Particles* _particles, 
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
            int _nLines,Boundary * _boundaryVector,
            Surfaces * _surfaces,
            Flags * _flags);

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const;
    
};
#endif