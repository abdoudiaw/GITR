#ifndef _REFLECTIVEBOUNDARYCONDITIONS_
#define _REFLECTIVEBOUNDARYCONDITIONS_

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

struct reflectiveBoundaryConditions {
    Particles * particles;
    int nLines;
    int nSurfaces;
    Boundary * boundaryVector;
    Flags * flags;

    reflectiveBoundaryConditions(Particles* _particles, 
            int _nLines,
            int _nSurfaces,
            Boundary * _boundaryVector,
            Flags * _flags);

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const;
    
};
#endif
