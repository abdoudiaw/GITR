#ifndef _GEOM_
#define _GEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <cmath>

#include "Boundary.h"
#include "Particles.h"
#include "surfaceReactions.h"
#include "pusher.h"


#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

__host__ __device__

struct geometryCheck {
  Particles *particlesPointer;
  const int nLines;
  Boundary *boundaryVector;
  Flags *flags;

  geometryCheck(Particles* _particlesPointer, int _nLines, const Boundary* const _boundaryVector, Flags* _flags);

  // Rest of the code remains the same

  __host__  __device__
  void operator()(std::size_t indx) const;
  void perform3DTetGeomCheck(std::size_t indx) const;
  // void updateParticlePosition(std::size_t indx, double x, double y, double z, bool flag) const;
  bool isIntersected(const gitr_precision p0[3], const gitr_precision p1[3], const Boundary& boundary, gitr_precision intersectionPoint[3]) const;
  bool isPointInsideTriangle( gitr_precision p[3], const Boundary& boundary) const;
  // void reflectDirection(const gitr_precision p1[3], const gitr_precision intersectionPoint[3], const gitr_precision normal[3], gitr_precision reflection_dir[3]) const; 
  void reflectDirection(const gitr_precision velocity[3], const gitr_precision intersectionPoint[3], const gitr_precision normal[3], gitr_precision reflection_dir[3]) const ;
  void reflectParticle(std::size_t indx, const gitr_precision reflection_point[3], const gitr_precision reflection_dir[3]) const;
void periodicBoundaryConditions(std::size_t indx,  gitr_precision temp_position_xyz[3]) const;
  void vectorSubtract(const gitr_precision a[3], const gitr_precision b[3], gitr_precision result[3]) const;
  gitr_precision vectorDotProduct(const gitr_precision a[3], const gitr_precision b[3]) const;
  void vectorScale(const gitr_precision a[3], const gitr_precision scalar, gitr_precision result[3]) const;
  void vectorNormalize(gitr_precision a[3]) const;
  gitr_precision vectorNorm(const gitr_precision a[3]) const;
  void vectorCrossProduct(const gitr_precision a[3], const gitr_precision b[3], gitr_precision result[3]) const;

void reflect3(gitr_precision v[3], const gitr_precision n[3]) const;
gitr_precision dot3(const gitr_precision v1[3], const gitr_precision v2[3]) const;

};

#endif
