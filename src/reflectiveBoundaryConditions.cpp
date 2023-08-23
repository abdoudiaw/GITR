#include <iostream>
#include <tuple>
#include <vector>
#include <libconfig.h++>
#include <netcdf>
#include "Boundary.h"
#include "reflectiveBoundaryConditions.h"


reflectiveBoundaryConditions::reflectiveBoundaryConditions(Particles* _particles,
 int _nLines, int _nSurfaces,
 Boundary* _boundaryVector, Flags* _flags) :
    particles(_particles), nLines(_nLines), nSurfaces(_nSurfaces), boundaryVector(_boundaryVector),  flags(_flags) { }

CUDA_CALLABLE_MEMBER_DEVICE
void reflectiveBoundaryConditions::operator()(std::size_t indx) const {

        gitr_precision thetaImpact = 0.0;
        gitr_precision particleTrackVector[3] = {0.0};
        gitr_precision surfaceNormalVector[3] = {0.0};
        gitr_precision vSampled[3] = {0.0};

        int wallHit = particles->surfaceHit[indx];
        int surfaceHit = boundaryVector[wallHit].surfaceNumber;
        int surface = boundaryVector[wallHit].surface;

        gitr_precision vx = particles->vx[indx];
        gitr_precision vy = particles->vy[indx];
        gitr_precision vz = particles->vz[indx];
        gitr_precision SURFACE_BUFFER = 1.0e-4;

        particleTrackVector[0] = vx;
        particleTrackVector[1] = vy;
        particleTrackVector[2] = vz;
        gitr_precision velocityNorm = std::sqrt(vx * vx + vy * vy + vz * vz);


        boundaryVector[wallHit].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);
        // printf(" surfaceNormalVector %g %g %g \n", surfaceNormalVector[0], surfaceNormalVector[1], surfaceNormalVector[2]);

                    
        particleTrackVector[0] = vx / velocityNorm;
        particleTrackVector[1] = vy / velocityNorm;
        particleTrackVector[2] = vz / velocityNorm;
        gitr_precision partDotNormal = particleTrackVector[0] * surfaceNormalVector[0] + particleTrackVector[1] * surfaceNormalVector[1] + particleTrackVector[2] * surfaceNormalVector[2];
        thetaImpact = std::acos(partDotNormal);

        vSampled[0] = vx;
        vSampled[1] = vy;
        vSampled[2] = vz;

        particles->vx[indx] = - vSampled[0];
        particles->vy[indx] = - vSampled[1];
        particles->vz[indx] = - vSampled[2];
        // print("surfaceNormalVector %g %g %g \n", surfaceNormalVector[0], surfaceNormalVector[1], surfaceNormalVector[2]);
        // exit(0);
        printf(" velocities %g %g %g \n", particles->vx[indx], particles->vy[indx], particles->vz[indx]);

        particles->xprevious[indx] = particles->x[indx] + surfaceNormalVector[0] * SURFACE_BUFFER;
        particles->yprevious[indx] = particles->y[indx] + surfaceNormalVector[1] * SURFACE_BUFFER;
        particles->zprevious[indx] = particles->z[indx] + surfaceNormalVector[2] * SURFACE_BUFFER;

        // print particle position
        printf(" particle position %g %g %g \n", particles->x[indx], particles->y[indx], particles->z[indx]);



        // // set periodic boundary conditions
        // gitr_precision xmax = 0.0919991;
        // gitr_precision xmin = 0;
        // gitr_precision ymax = 0.0919991;
        // gitr_precision ymin = 0;
        // gitr_precision zmax = 0.0919991;
        // gitr_precision zmin = 0;

        // if (particles->x[indx] > xmax) {
        //     particles->x[indx] = xmin + (particles->x[indx] - xmax);
        // }
        // else if (particles->x[indx] < xmin) {
        //     particles->x[indx] = xmax - (xmin - particles->x[indx]);
        // }

        // if (particles->y[indx] > ymax) {
        //     particles->y[indx] = ymin + (particles->y[indx] - ymax);
        // }
        // else if (particles->y[indx] < ymin) {
        //     particles->y[indx] = ymax - (ymin - particles->y[indx]);
        // }

        // if (particles->z[indx] > zmax) {
        //     particles->z[indx] = zmin + (particles->z[indx] - zmax);
        // }
        // else if (particles->z[indx] < zmin) {
        //     particles->z[indx] = zmax - (zmin - particles->z[indx]);
        // }

}

 