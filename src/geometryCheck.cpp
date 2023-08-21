#include "geometryCheck.h"
#include "flags.h"

geometryCheck::geometryCheck(Particles* _particlesPointer, int _nLines, const Boundary* const _boundaryVector, Flags* _flags)
    : particlesPointer(_particlesPointer), nLines(_nLines), boundaryVector(const_cast<Boundary*>(_boundaryVector)), flags(_flags)
{}

__host__ __device__
void geometryCheck::operator()(std::size_t indx) const {
  if (particlesPointer->hitWall[indx] == 0.0) {
    perform3DTetGeomCheck(indx);
  }
}

__host__ __device__
void geometryCheck::perform3DTetGeomCheck(std::size_t indx) const {
  if (flags->USE3DTETGEOM <= 0) {
    return;
  }

  int hitSurface = 0;
  gitr_precision x = particlesPointer->x[indx], y = particlesPointer->y[indx], z = particlesPointer->z[indx];
  gitr_precision xprev = particlesPointer->xprevious[indx], yprev = particlesPointer->yprevious[indx], zprev = particlesPointer->zprevious[indx];
  gitr_precision dpath = sqrt((x - xprev) * (x - xprev) + (y - yprev) * (y - yprev) + (z - zprev) * (z - zprev));
  gitr_precision p0[3] = {particlesPointer->xprevious[indx], particlesPointer->yprevious[indx], particlesPointer->zprevious[indx]};
  gitr_precision p1[3] = {particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx]};
  int nBoundariesCrossed = 0, boundariesCrossed[6] = {0};
  gitr_precision nearest_boundary_distance = 1.0e20, temp_position_xyz[3] = {0.0};
  int nearest_boundary_index = 0;

  for (int i = 0; i < nLines; i++) {
    gitr_precision intersectionPoint[3] = {0.0};
    if (isIntersected(p0, p1, boundaryVector[i], intersectionPoint)) {
      if (isPointInsideTriangle(intersectionPoint, boundaryVector[i])) {

        gitr_precision dist_to_p = sqrt((intersectionPoint[0] - xprev) * (intersectionPoint[0] - xprev) + (intersectionPoint[1] - yprev) * (intersectionPoint[1] - yprev) + (intersectionPoint[2] - zprev) * (intersectionPoint[2] - zprev));
        if (dist_to_p < nearest_boundary_distance) {
          nearest_boundary_distance = dist_to_p;
          nBoundariesCrossed++;
          nearest_boundary_index = i;
          for (int j = 0; j < 3; j++) {
            temp_position_xyz[j] = intersectionPoint[j];
          }
        }
      }
      gitr_precision normal[3] = {
          boundaryVector[static_cast<size_t>(nearest_boundary_index)].a,
          boundaryVector[static_cast<size_t>(nearest_boundary_index)].b,
          boundaryVector[static_cast<size_t>(nearest_boundary_index)].c
      };
    }
  }

  if (nBoundariesCrossed == 0) {
    particlesPointer->distTraveled[indx] += dpath;
  } else {
    particlesPointer->hitWall[indx] = 1.0;
    gitr_precision normal[3] = {boundaryVector[nearest_boundary_index].a, boundaryVector[nearest_boundary_index].b, boundaryVector[nearest_boundary_index].c};
    gitr_precision reflection_dir[3];
    gitr_precision velocity[3];
    // velocity[0] = particlesPointer->vx[indx];
    // velocity[1] = particlesPointer->vy[indx];
    // velocity[2] = particlesPointer->vz[indx];
    // reflectDirection(velocity, temp_position_xyz, normal, reflection_dir);
    // reflectParticle(indx, temp_position_xyz, reflection_dir);

    // call periodic boundary conditions
    // periodicBoundaryConditions(indx, temp_position_xyz);
                       // apply periodic boundary conditions
                       int nP = 20;
                    gitr_precision eps = 1e-6;
        gitr_precision xmin = 0;
        gitr_precision xmax = 0.0919991;
        gitr_precision ymin = 0;
        gitr_precision ymax = 0.0919991;
        gitr_precision zmin = 0;
        gitr_precision zmax = 0.0919991;

              for (int part=0; part < nP; part++) {
                  if (particlesPointer->x[part] < xmin) {
                      particlesPointer->x[part] = xmax - eps;
                  }
                  if (particlesPointer->x[part] > xmax) {
                      particlesPointer->x[part] = xmin + eps;
                  }
                  if (particlesPointer->y[part] < ymin) {
                      particlesPointer->y[part] = ymax - eps;
                  }
                  if (particlesPointer->y[part] > ymax) {
                      particlesPointer->y[part] = ymin + eps;
                  }
                  if (particlesPointer->z[part] < zmin) {
                      particlesPointer->z[part] = zmax - eps;
                  }
                  if (particlesPointer->z[part] > zmax) {
                      particlesPointer->z[part] = zmin + eps;
                  }

              };


    particlesPointer->distTraveled[indx] += nearest_boundary_distance;
    particlesPointer->surfaceHit[indx] = nearest_boundary_index;
  }
}

__host__ __device__
void geometryCheck::reflectDirection(const gitr_precision velocity[3], const gitr_precision intersectionPoint[3], const gitr_precision normal[3], gitr_precision reflection_dir[3]) const {
  gitr_precision incident_dir[3];
  gitr_precision dot = vectorDotProduct(velocity, normal);
  vectorScale(normal, 2.0 * dot, reflection_dir);
  vectorSubtract(velocity, reflection_dir, reflection_dir);
  vectorNormalize(reflection_dir);
}

__host__ __device__
void geometryCheck::reflectParticle(std::size_t indx, const gitr_precision reflection_point[3], const gitr_precision reflection_dir[3]) const {
    gitr_precision small_step = 1e-06;  // Small step to move the particle away from the wall
    particlesPointer->x[indx] = reflection_point[0]; //  - small_step * reflection_dir[0];
    particlesPointer->y[indx] = reflection_point[1]; // - small_step * reflection_dir[1];
    particlesPointer->z[indx] = reflection_point[2]; // - small_step * reflection_dir[2];

    gitr_precision particle_velocity[3] = {particlesPointer->vx[indx], particlesPointer->vy[indx], particlesPointer->vz[indx]};
    gitr_precision particle_velocity_magnitude = vectorNorm(particle_velocity);

    // refelct particle usign reflect3
    reflect3(particle_velocity, reflection_dir);
    // update particle velocities
    particlesPointer->vx[indx] = particle_velocity[0];
    particlesPointer->vy[indx] = particle_velocity[1];
    particlesPointer->vz[indx] = particle_velocity[2];
  }

__host__ __device__
void geometryCheck::reflect3(gitr_precision v[3], const gitr_precision n[3]) const
{
  double dot = dot3(v,n);
  v[0] -= 2.0*dot*n[0];
  v[1] -= 2.0*dot*n[1];
  v[2] -= 2.0*dot*n[2];
}

void geometryCheck::periodicBoundaryConditions(std::size_t indx,  gitr_precision temp_position_xyz[3]) const {
        gitr_precision eps = 1e-6;
        gitr_precision xmin = 0;
        gitr_precision xmax = 0.0919991;
        gitr_precision ymin = 0;
        gitr_precision ymax = 0.0919991;
        gitr_precision zmin = 0;
        gitr_precision zmax = 0.0919991;

        if (temp_position_xyz[0] < xmin) {
          temp_position_xyz[0] = xmax - eps;
        }
        if (temp_position_xyz[0] > xmax) {
          temp_position_xyz[0] = xmin + eps;
        }
        if (temp_position_xyz[1] < ymin) {
          temp_position_xyz[1] = ymax - eps;
        }
        if (temp_position_xyz[1] > ymax) {
          temp_position_xyz[1] = ymin + eps;
        }
        if (temp_position_xyz[2] < zmin) {
          temp_position_xyz[2] = zmax - eps;
        }
        if (temp_position_xyz[2] > zmax) {
          temp_position_xyz[2] = zmin + eps;
        }
        particlesPointer->x[indx] = temp_position_xyz[0];
        particlesPointer->y[indx] = temp_position_xyz[1];
        particlesPointer->z[indx] = temp_position_xyz[2];

}
  

__host__ __device__
gitr_precision geometryCheck::dot3(const gitr_precision v1[3], const gitr_precision v2[3]) const
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


__host__ __device__
bool geometryCheck::isIntersected(const gitr_precision p0[3], const gitr_precision p1[3], const Boundary& boundary, gitr_precision intersectionPoint[3]) const {
  gitr_precision a = boundary.a;
  gitr_precision b = boundary.b;
  gitr_precision c = boundary.c;
  gitr_precision d = boundary.d;
  gitr_precision plane_norm = boundary.plane_norm;

  gitr_precision pointToPlaneDistance0 = (a * p0[0] + b * p0[1] + c * p0[2] + d) / plane_norm;
  gitr_precision pointToPlaneDistance1 = (a * p1[0] + b * p1[1] + c * p1[2] + d) / plane_norm;

  gitr_precision signPoint0 = copysign(1.0, pointToPlaneDistance0);
  gitr_precision signPoint1 = copysign(1.0, pointToPlaneDistance1);

  if (signPoint0 != signPoint1) {
    gitr_precision t = -(a * p0[0] + b * p0[1] + c * p0[2] + d) / (a * (p1[0] - p0[0]) + b * (p1[1] - p0[1]) + c * (p1[2] - p0[2]));

    for (int i = 0; i < 3; i++) {
      intersectionPoint[i] = p0[i] + t * (p1[i] - p0[i]);
    }
    // printf(" intersection true!");
    return true;
  }
  // printf(" intersection false!");
  return false;
}

__host__ __device__
bool geometryCheck::isPointInsideTriangle(gitr_precision p[3], const Boundary& boundary) const {
  gitr_precision A[3] = {boundary.x1, boundary.y1, boundary.z1};
  gitr_precision B[3] = {boundary.x2, boundary.y2, boundary.z2};
  gitr_precision C[3] = {boundary.x3, boundary.y3, boundary.z3};

  gitr_precision AB[3], AC[3], BC[3], CA[3];
  gitr_precision Ap[3], Bp[3], Cp[3];
  gitr_precision crossABAp[3], crossBCBp[3], crossCACp[3];
  gitr_precision normalVector[3];

  vectorSubtract(B, A, AB);
  vectorSubtract(C, A, AC);
  vectorSubtract(C, B, BC);
  vectorSubtract(A, C, CA);

  vectorSubtract(p, A, Ap);
  vectorSubtract(p, B, Bp);
  vectorSubtract(p, C, Cp);

  vectorCrossProduct(AB, AC, normalVector);
  vectorCrossProduct(AB, Ap, crossABAp);
  vectorCrossProduct(BC, Bp, crossBCBp);
  vectorCrossProduct(CA, Cp, crossCACp);

  gitr_precision signDot0 = copysign(1.0, vectorDotProduct(crossABAp, normalVector));
  gitr_precision signDot1 = copysign(1.0, vectorDotProduct(crossBCBp, normalVector));
  gitr_precision signDot2 = copysign(1.0, vectorDotProduct(crossCACp, normalVector));
  gitr_precision totalSigns = 1.0 * abs(signDot0 + signDot1 + signDot2);

  if (totalSigns == 3.0 && vectorNorm(crossABAp) != 0.0 && vectorNorm(crossBCBp) != 0.0 && vectorNorm(crossCACp) != 0.0) {
    return true;
  }

  return false;
}

__host__ __device__
void geometryCheck::vectorSubtract(const gitr_precision a[3], const gitr_precision b[3], gitr_precision result[3]) const {
  for (int i = 0; i < 3; i++) {
    result[i] = a[i] - b[i];
  }
}

__host__ __device__
void geometryCheck::vectorScale(const gitr_precision a[3], const gitr_precision scalar, gitr_precision result[3]) const {
  for (int i = 0; i < 3; i++) {
    result[i] = a[i] * scalar;
  }
}

__host__ __device__
void geometryCheck::vectorNormalize(gitr_precision a[3]) const {
  gitr_precision norm = vectorNorm(a);
  if (norm != 0.0) {
    for (int i = 0; i < 3; i++) {
      a[i] /= norm;
    }
  }
}

__host__ __device__
gitr_precision geometryCheck::vectorNorm(const gitr_precision a[3]) const {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

__host__ __device__
gitr_precision geometryCheck::vectorDotProduct(const gitr_precision a[3], const gitr_precision b[3]) const {
  gitr_precision dot = 0.0;
  for (int i = 0; i < 3; i++) {
    dot += a[i] * b[i];
  }
  return dot;
}

__host__ __device__
void geometryCheck::vectorCrossProduct(const gitr_precision a[3], const gitr_precision b[3], gitr_precision result[3]) const {
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];
}
