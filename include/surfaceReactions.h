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
                         const std::vector<double>& gridx, const std::vector<double>& gridz, const std::vector<double>& data);
struct surfaceReactions {
    Particles * particles;
    int nLines;
    int nSurfaces;
    Boundary * boundaryVector;
    Flags * flags;

    surfaceReactions(Particles* _particles, 
            int _nLines, int _nSurfaces,
            Boundary * _boundaryVector,
            Flags * _flags);

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const;
    
};
#endif
