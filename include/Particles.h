#ifndef _PARTICLES_
#define _PARTICLES_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#include "flags.h"
#include <libconfig.h++>

#include <cstdlib>
#include <stdio.h>

#ifdef __CUDACC__
#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#endif

#include <random>
#include <math.h>
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

class Particles : public ManagedAllocation {
public:
  std::size_t nParticles;
  sim::Array<int> index;
  sim::Array<gitr_precision> x;
  sim::Array<gitr_precision> y;
  sim::Array<gitr_precision> z;
  sim::Array<gitr_precision> xprevious;
  sim::Array<gitr_precision> yprevious;
  sim::Array<gitr_precision> zprevious;
  sim::Array<gitr_precision> v;
  sim::Array<gitr_precision> vx;
  sim::Array<gitr_precision> vy;
  sim::Array<gitr_precision> vz;
  sim::Array<gitr_precision> Z;
  sim::Array<gitr_precision> amu;
  sim::Array<gitr_precision> charge;
  sim::Array<gitr_precision> newVelocity;
  sim::Array<int> tt;
  sim::Array<int> hasLeaked;
  sim::Array<gitr_precision> leakZ;
#ifdef __CUDACC__
  sim::Array<curandState> stream_ionize;
#else
  sim::Array<std::mt19937> stream_ionize;
#endif

  sim::Array<gitr_precision> hitWall;
  sim::Array<int> surfaceHit;
  sim::Array<int> firstCollision;
  sim::Array<gitr_precision> transitTime;
  sim::Array<gitr_precision> distTraveled;
  sim::Array<int> wallIndex;
  sim::Array<gitr_precision> distanceTraveled;
  sim::Array<gitr_precision> weight;
  sim::Array<gitr_precision> dt;
  sim::Array<gitr_precision> time;
  sim::Array<bool> advance;
  std::vector<std::string> materialName; 

  CUDA_CALLABLE_MEMBER
  void setParticle(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                   gitr_precision Ex, gitr_precision Ey, gitr_precision Ez,
                   gitr_precision Z, gitr_precision amu, gitr_precision charge,
                   const std::string& materialName)
  {

    // this->index[indx] = indx;
    this->xprevious[indx] = x;
    this->yprevious[indx] = y;
    this->zprevious[indx] = z;
    this->x[indx] = x;
    this->y[indx] = y;
    this->z[indx] = z;
    this->Z[indx] = Z;
    this->charge[indx] = charge;
    this->amu[indx] = amu;
    this->hitWall[indx] = 0.0;
    this->wallIndex[indx] = 0;
    this->vx[indx] = 0.0;
    this->vy[indx] = 0.0;
    this->vz[indx] = 0.0;
    this->v[indx] = 0.0;
  };

  CUDA_CALLABLE_MEMBER
  void setParticleV(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                    gitr_precision Vx, gitr_precision Vy, gitr_precision Vz,
                    gitr_precision Z, gitr_precision amu, gitr_precision charge,
                    gitr_precision dt, const std::string& materialName)
  {
    int indTmp = indx;
    this->index[indx] = indTmp;
    this->xprevious[indx] = x;
    this->yprevious[indx] = y;
    this->zprevious[indx] = z;
    this->x[indx] = x;
    this->y[indx] = y;
    this->z[indx] = z;
    this->Z[indx] = Z;
    this->charge[indx] = charge;
    this->amu[indx] = amu;
    this->hitWall[indx] = 0.0;
    this->wallIndex[indx] = 0;
    this->vx[indx] = Vx;
    this->vy[indx] = Vy;
    this->vz[indx] = Vz;
    this->v[indx] = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
    this->dt[indx] = dt;
    this->materialName[indx] = materialName;

  };


  CUDA_CALLABLE_MEMBER
  Particles(std::size_t nP,std::size_t nStreams,libconfig::Config &cfg, Flags *gitr_flags) :
    nParticles{nP},
    index{nParticles, 0},
    x{nP,0.0},
    y{nP,0.0},
    z{nP,0.0},
    xprevious{nParticles,0.0},
    yprevious{nParticles,0.0},
    zprevious{nParticles,0.0},
    v{nParticles, 0.0},
    vx{nParticles,0.0},
    vy{nParticles,0.0},
    vz{nParticles,0.0},
    Z{nParticles,0.0},
    amu{nParticles,0.0},
    charge{nParticles,0.0},
    newVelocity{nParticles},
    tt{nParticles, 0},
    hasLeaked{nParticles, 0},
    leakZ{nParticles,0.0},
    stream_ionize{nParticles},
    hitWall{nParticles, 0.0},
    surfaceHit{nParticles, -1},
    firstCollision{nParticles, 1},
    transitTime{nParticles, 0.0},
    distTraveled{nParticles, 0.0},
    wallIndex{nParticles,0},
    distanceTraveled{nParticles,0.0},
    weight{nParticles, 1.0},
    dt{nParticles,0.0},
    time{nParticles,0.0},
    materialName{nParticles, "None"},
    advance{nParticles,false} {};
};



#endif
