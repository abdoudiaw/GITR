#ifndef _BOUNDARYCONDITIONS_
#define _BOUNDARYCONDITIONS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
#include <stdio.h>
#include "array.h"

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#else
#include <random>
#endif

#include "Particles.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

class boundaryConditions {
public:
    boundaryConditions(Particles *_particlesPointer);
    CUDA_CALLABLE_MEMBER void operator()(std::size_t indx);

private:
    Particles *particlesPointer;
};

boundaryConditions::boundaryConditions(Particles *_particlesPointer)
    : particlesPointer(_particlesPointer)
{}

CUDA_CALLABLE_MEMBER    
void boundaryConditions::operator()(std::size_t indx) {
    gitr_precision position[3] = {particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx]};

    gitr_precision eps = 1e-6;
    gitr_precision xmin = 0;
    gitr_precision xmax = 0.0919991;
    gitr_precision ymin = 0;
    gitr_precision ymax = 0.0919991;
    gitr_precision zmin = 0;
    gitr_precision zmax = 0.0919991;

    // apply conditions
    if (position[0] < xmin) {
        position[0] = xmax - eps;
    }
    if (position[0] > xmax) {
        position[0] = xmin + eps;
    }
    if (position[1] < ymin) {
        position[1] = ymax - eps;
    }
    if (position[1] > ymax) {
        position[1] = ymin + eps;
    }
    if (position[2] < zmin) {
        position[2] = zmax - eps;
    }
    if (position[2] > zmax) {
        position[2] = zmin + eps;
    }

    // update position
    particlesPointer->x[indx] = position[0];
    particlesPointer->y[indx] = position[1];
    particlesPointer->z[indx] = position[2];
}

#endif
