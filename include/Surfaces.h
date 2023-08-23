#ifndef _SURFACES_
#define _SURFACES_

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
#include <random>
#include <string>
#include "array.h"

#if USE_DOUBLE
using gitr_precision = double;
#else
using gitr_precision = float;
#endif

class Surfaces : public ManagedAllocation {
public:
    int nSurfaces;  
    int nE;
    int nA;


    sim::Array<gitr_precision> sumParticlesStrike;
    sim::Array<gitr_precision> sumWeightStrike;
    sim::Array<gitr_precision> grossDeposition;
    sim::Array<gitr_precision> grossErosion;
    sim::Array<gitr_precision> aveSputtYld;
    sim::Array<gitr_precision> sputtYldCount;
    sim::Array<gitr_precision> energyDistribution;
    sim::Array<gitr_precision> sputtDistribution;
    sim::Array<gitr_precision> reflDistribution;
    sim::Array<gitr_precision> gridE;
    sim::Array<gitr_precision> gridA;

    std::string materialName;

    CUDA_CALLABLE_MEMBER
    Surfaces(std::size_t nSurfaces, std::size_t nE, std::size_t nA) :
        sumParticlesStrike{nSurfaces, 0}, gridE{nE, 0.0}, gridA{nA, 0.0},
        sumWeightStrike{nSurfaces, 0.0}, grossDeposition{nSurfaces, 0.0},
        grossErosion{nSurfaces, 0.0}, aveSputtYld{nSurfaces, 0.0}, sputtYldCount{nSurfaces, 0},
        energyDistribution{nSurfaces * nE * nA, 0.0}, sputtDistribution{nSurfaces * nE * nA, 0.0},
        reflDistribution{nSurfaces * nE * nA, 0.0} {};
};

#endif
