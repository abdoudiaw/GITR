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
#include "boundaryConditions.h"


#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

class boundaryConditions {
public:
    boundaryConditions(Particles *_particlesPointer, const Domain& domain);
    CUDA_CALLABLE_MEMBER void operator()(std::size_t indx);

private:
    Particles *particlesPointer;
    Domain domain;
};

boundaryConditions::boundaryConditions(Particles *_particlesPointer, const Domain& _domain)
    : particlesPointer(_particlesPointer), domain(_domain)
{}

CUDA_CALLABLE_MEMBER    
void boundaryConditions::operator()(std::size_t indx) {
    gitr_precision position[3] = {particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx]};

    gitr_precision eps = 1e-6;
    gitr_precision xmin =domain.xmin;
    gitr_precision xmax =domain.xmax;
    gitr_precision ymin =domain.ymin;
    gitr_precision ymax =domain.ymax;
    gitr_precision zmin =domain.zmin;
    gitr_precision zmax =domain.zmax;
    

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