
#include "Particles.h"


#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#include <curand_kernel.h>
#include <thrust/random.h>
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#ifdef __GNUC__
#include <random>
#endif

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

#if USE_CUDA
CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(curandState *state,int indx);
double get_rand_double(curandState *state,int indx);
#else
CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(std::mt19937 *state,int indx);
double get_rand_double(std::mt19937 *state,int indx);
#endif


std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>> process_rates(int charge);
gitr_precision rateCoeffInterp(int charge, gitr_precision te, gitr_precision ne, int nT, int nD, const std::vector<double>& rateGrid_Tempp, const std::vector<double>& rateGrid_Densp, const std::vector<double>& Ratesp);

template <typename T=std::mt19937>
struct elementaryProcesses {
  Particles *particlesPointer;
  int nR_Dens;
  int nZ_Dens;
  gitr_precision *DensGridr;
  gitr_precision *DensGridz;
  gitr_precision *ne;
  int nR_Temp;
  int nZ_Temp;
  gitr_precision *TempGridr;
  gitr_precision *TempGridz;
  gitr_precision *te;
  gitr_precision tion;
  void (elementaryProcesses::*func)(std::size_t);
  int xx1;
  T *state;
  Flags *flags;
  
  elementaryProcesses( Particles *_particlesPointer, 
  T *_state,
         int _nR_Dens, int _nZ_Dens, gitr_precision *_DensGridr, gitr_precision *_DensGridz,
         gitr_precision *_ne, int _nR_Temp, int _nZ_Temp, gitr_precision *_TempGridr,
         gitr_precision *_TempGridz, gitr_precision *_te, Flags *_flags);

  CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx);
};

// END ELEMENTARY_PROCESSES_H