#ifndef _FLAGS_
#define _FLAGS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "libconfig.h++"
#include "utils.h"

class Flags : public ManagedAllocation
{
  private:
  public:
    const bool USE3DTETGEOM;
    const bool USE_SURFACEMODEL;
    const bool USE_SHEATHEFIELD;
        const bool USE_IONIZATION;
    const bool USE_RECOMBINATION;
   CUDA_CALLABLE_MEMBER
   Flags(libconfig::Config &cfg) : 
    USE3DTETGEOM{initialize_flags(cfg,"USE3DTETGEOM")},
    USE_SURFACEMODEL{initialize_flags(cfg,"USE_SURFACEMODEL")},
        USE_IONIZATION{initialize_flags(cfg,"USE_IONIZATION")},
    USE_RECOMBINATION{initialize_flags(cfg,"USE_RECOMBINATION")},
    USE_SHEATHEFIELD{initialize_flags(cfg,"USE_SHEATHEFIELD")} {};
   bool initialize_flags(libconfig::Config &cfg, std::string s);
};
#endif
