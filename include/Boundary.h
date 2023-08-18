#ifndef _BOUNDARY_
#define _BOUNDARY_

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

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

class Boundary 
{
  public:
    int periodic;
    int periodic_bc_x;
    int pointLine;
    int surfaceNumber;
    int surface;
    int inDir;
    gitr_precision x1;
    gitr_precision y1;
    gitr_precision z1;
    gitr_precision x2;
    gitr_precision y2;
    gitr_precision z2;
    gitr_precision a;
    gitr_precision b;
    gitr_precision c;
    gitr_precision d;
    gitr_precision plane_norm; //16
    gitr_precision x3;
    gitr_precision y3;
    gitr_precision z3;
    gitr_precision area;
    gitr_precision slope_dzdx;
    gitr_precision intercept_z;
    gitr_precision periodic_bc_x0;    
    gitr_precision periodic_bc_x1;    
    gitr_precision Z;
    gitr_precision amu;
    gitr_precision potential;
    gitr_precision ChildLangmuirDist;

    gitr_precision hitWall;
    gitr_precision length;
    gitr_precision distanceToParticle;
    gitr_precision angle;
    gitr_precision fd;
    gitr_precision density;
    gitr_precision ti;
    gitr_precision ne;
    gitr_precision te;
    gitr_precision debyeLength;
    gitr_precision larmorRadius;
    gitr_precision flux;
    gitr_precision startingParticles;
    gitr_precision impacts;
    gitr_precision redeposit;
    gitr_precision unit_vec0;
    gitr_precision unit_vec1;
    gitr_precision unit_vec2;

    CUDA_CALLABLE_MEMBER
    void getSurfaceParallel(gitr_precision A[],gitr_precision y,gitr_precision x,
                            int use_3d_geom )
    {
    gitr_precision norm;
    if( use_3d_geom > 0 )
    {
    norm = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    A[1] = (y2 - y1) / norm;
    }
    else
    {
    norm = std::sqrt((x2 - x1) * (x2 - x1) + (z2 - z1) * (z2 - z1));
    A[1] = 0.0;
    }
        A[0] = (x2-x1)/norm;
        A[2] = (z2-z1)/norm;

    }

  CUDA_CALLABLE_MEMBER
  void getSurfaceNormal(gitr_precision B[], gitr_precision y, gitr_precision x, int use_3d_geom ) {
    if( use_3d_geom > 0 )
    {
    B[0] = a / plane_norm;
    B[1] = b / plane_norm;
    B[2] = c / plane_norm;
    }
    else
    {
    gitr_precision perpSlope = 0.0;
    if (slope_dzdx == 0.0) {
      perpSlope = 1.0e12;
    } else {
      perpSlope = -std::copysign(1.0, slope_dzdx) / std::abs(slope_dzdx);
    }
    gitr_precision Br = 1.0 / std::sqrt(perpSlope * perpSlope + 1.0);
    gitr_precision Bt = 0.0;
    B[2] = std::copysign(1.0,perpSlope) * std::sqrt(1 - Br * Br);
      B[0] = Br;
      B[1] = Bt;
    
    }
    }
    CUDA_CALLABLE_MEMBER
        void transformToSurface(gitr_precision C[],gitr_precision y, gitr_precision x, int use_3d_geom)
        {
            gitr_precision X[3] = {0.0};
            gitr_precision Y[3] = {0.0};
            gitr_precision Z[3] = {0.0};
            gitr_precision tmp[3] = {0.0};
            
            getSurfaceParallel(X,y,x, use_3d_geom );
            getSurfaceNormal(Z,y,x, use_3d_geom );

            Y[0] = Z[1]*X[2] - Z[2]*X[1]; 
            Y[1] = Z[2]*X[0] - Z[0]*X[2]; 
            Y[2] = Z[0]*X[1] - Z[1]*X[0];

            tmp[0] = X[0]*C[0] + Y[0]*C[1] + Z[0]*C[2];
            tmp[1] = X[1]*C[0] + Y[1]*C[1] + Z[1]*C[2];
            tmp[2] = X[2]*C[0] + Y[2]*C[1] + Z[2]*C[2];

            C[0] = tmp[0];
            C[1] = tmp[1];
            C[2] = tmp[2];

        }
};
#endif
