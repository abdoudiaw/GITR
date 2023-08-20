#ifndef _BORIS_
#define _BORIS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#include "thrust/extrema.h"
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <algorithm>
#include "Particles.h"
#include "Boundary.h"
#include "interpolater.h"
#include "flags.h"
#include <algorithm>

#include <tuple>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

CUDA_CALLABLE_MEMBER
void vectorAdd(gitr_precision A[], gitr_precision B[],gitr_precision C[]);

CUDA_CALLABLE_MEMBER
void vectorSubtract(gitr_precision A[], gitr_precision B[],gitr_precision C[]);

CUDA_CALLABLE_MEMBER
void vectorScalarMult(gitr_precision a, gitr_precision B[],gitr_precision C[]);

CUDA_CALLABLE_MEMBER
void vectorAssign(gitr_precision a, gitr_precision b,gitr_precision c, gitr_precision D[]);

CUDA_CALLABLE_MEMBER
gitr_precision vectorNorm(gitr_precision A[]);

CUDA_CALLABLE_MEMBER
void vectorNormalize(gitr_precision A[],gitr_precision B[]);

CUDA_CALLABLE_MEMBER
gitr_precision vectorDotProduct(gitr_precision A[], gitr_precision B[]);

CUDA_CALLABLE_MEMBER
void vectorCrossProduct(gitr_precision A[], gitr_precision B[], gitr_precision C[]);

CUDA_CALLABLE_MEMBER
gitr_precision getSheathElectricField ( gitr_precision x0, gitr_precision y, gitr_precision z, gitr_precision E[], Boundary *boundaryVector, int nLines,
         int&  closestBoundaryIndex, int use_3d_geom, gitr_precision mass, gitr_precision charge, gitr_precision Bmagnitude);
CUDA_CALLABLE_MEMBER
gitr_precision compute_cross_product(const gitr_precision A[3], const gitr_precision B[3], int i); 

// Function declaration for boris_push
CUDA_CALLABLE_MEMBER
std::tuple<gitr_precision, gitr_precision, gitr_precision,
           gitr_precision, gitr_precision, gitr_precision> 
borisPush(Particles *particlesPointer, std::size_t indx,
          gitr_precision E[3], gitr_precision B[3]);

struct pusher{ 
    Particles *particlesPointer;
    Boundary *boundaryVector;
    int nR_Bfield;
    int nZ_Bfield;
    gitr_precision * BfieldGridRDevicePointer;
    gitr_precision * BfieldGridZDevicePointer;
    gitr_precision * BfieldRDevicePointer;
    gitr_precision * BfieldZDevicePointer;
    gitr_precision * BfieldTDevicePointer;
    Flags* gitr_flags;
    int nLines;

    pusher(Particles *_particlesPointer , Boundary *_boundaryVector,int _nLines,
            int _nR_Bfield, int _nZ_Bfield,
            gitr_precision * _BfieldGridRDevicePointer,
            gitr_precision * _BfieldGridZDevicePointer,
            gitr_precision * _BfieldRDevicePointer,
            gitr_precision * _BfieldZDevicePointer,
            gitr_precision * _BfieldTDevicePointer,
            Flags* _gitr_flags);

    CUDA_CALLABLE_MEMBER    
    void operator()(std::size_t indx); 
};
#endif
