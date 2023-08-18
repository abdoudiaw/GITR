#ifndef _COLLISIONS_
#define _COLLISIONS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include <cmath>
#include "interpolater.h"
#include "pusherBoris.h"
#include "array.h"
#include "flags.h"
#include <cmath>
#include <limits>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#endif


// Compute the index of the cell based on particle coordinates
int computeCellIndex(gitr_precision x, gitr_precision y, gitr_precision z, const Boundary* boundaryVector, int nLines) {
    int cellIndex = 0;
    int numberOfCellsx = 20;
    int numberOfCellsy = 20;
    int numberOfCellsz = 20;

    // Compute the grid dimensions
    gitr_precision gridMinX = boundaryVector[0].x;
    gitr_precision gridMinY = boundaryVector[0].y;
    gitr_precision gridMinZ = boundaryVector[0].z;

    std::cout << "gridMinX: " << gridMinX << std::endl;
    std::cout << "gridMinY: " << gridMinY << std::endl;
    std::cout << "gridMinZ: " << gridMinZ << std::endl;

    // Compute the cell index based on the coordinates using a uniform grid
    int cellIndexX = static_cast<int>((x - gridMinX) / cellSizeX);
    int cellIndexY = static_cast<int>((y - gridMinY) / cellSizeY);
    int cellIndexZ = static_cast<int>((z - gridMinZ) / cellSizeZ);

    // Validate the cell indices to ensure they are within the grid dimensions
    cellIndexX = std::clamp(cellIndexX, 0, numberOfCellsx - 1);
    cellIndexY = std::clamp(cellIndexY, 0, numberOfCellsy - 1);
    cellIndexZ = std::clamp(cellIndexZ, 0, numberOfCellsz - 1);

    // Compute the cell index using a linear indexing scheme
    cellIndex = cellIndexX + cellIndexY * numberOfCellsx + cellIndexZ * (numberOfCellsx * numberOfCellsy);

    return cellIndex;
}


// Get a random particle index within the same cell
int getRandomParticleIndexInCell(int cellIndex) {
  int randomIndx = -1;
  // Get a random particle index within the specified cell
  // You need to define your own logic here based on your grid representation
  // ...
  return randomIndx;
}


gitr_precision  CoulombLogarithm(gitr_precision ne, gitr_precision ti_eV, gitr_precision Z1, gitr_precision Z2) {
  // note - this assumes ion and electron temperatures are the same

  if (ne == 0 || ti_eV == 0) {
    return 2;
  }

  double e2_cgs = 1.44e-7;  // eV-cm

  gitr_precision lam_e = 4 * M_PI * e2_cgs * ne * 1e-6 / ti_eV;  // 1/cm^2, this is really 1/lam_e
  gitr_precision g_12 = Z1 * Z2 * e2_cgs * sqrt(lam_e) / ti_eV;

  gitr_precision K11;
  if (g_12 < 1.) {
    K11 = -0.25 * log(1.466 * g_12 - 1.7836 * pow(g_12, 2.0) + 1.4313 * pow(g_12, 3.0) - 0.55833 * pow(g_12, 4) +
                      0.061162 * pow(g_12, 5)) *
          2.0;
  } else {
    K11 = (0.081033 - 0.091336 * log(g_12) + 0.05176 * pow(log(g_12), 2)) /
          (1.0 - 0.50026 * g_12 + 0.17044 * pow(g_12, 2)) * 2.0;
  }

  return fmax(K11, 2.0);
}


struct collisions {
  Particles* particlesPointer;
  gitr_precision dt;
  int Ncell; // Number of cells in the grid
  gitr_precision* cellDensity; // Array storing local densities for each cell
  gitr_precision* cellTemperature; // Array storing local temperatures for each cell
  gitr_precision background_Z;
  gitr_precision background_amu;

#ifdef __CUDACC__
  curandState* state;
#else
  std::mt19937* state;
#endif

  collisions(Particles* _particlesPointer, gitr_precision _dt,
#ifdef __CUDACC__
                   curandState* _state,
#else
                   std::mt19937* _state,
#endif
                   int _Ncell, gitr_precision* _cellDensity, gitr_precision* _cellTemperature,
                   gitr_precision _background_Z, gitr_precision _background_amu)
      : particlesPointer(_particlesPointer),
        dt(_dt),
        Ncell(_Ncell),
        cellDensity(_cellDensity),
        cellTemperature(_cellTemperature),
        background_Z(_background_Z),
        background_amu(_background_amu),
        state(_state) {}

  CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx) {
    if (particlesPointer->hitWall[indx] == 0.0 && particlesPointer->charge[indx] != 0.0) {
      gitr_precision x = particlesPointer->xprevious[indx];
      gitr_precision y = particlesPointer->yprevious[indx];
      gitr_precision z = particlesPointer->zprevious[indx];
      gitr_precision vx = particlesPointer->vx[indx];
      gitr_precision vy = particlesPointer->vy[indx];
      gitr_precision vz = particlesPointer->vz[indx];

      int cellIndex = computeCellIndex(x, y, z); // Compute the index of the cell for the particle

#ifdef __CUDACC__

      gitr_precision n1 = curand_normal(&state[indx]);
      gitr_precision n2 = curand_normal(&state[indx]);
      gitr_precision r1 = curand_uniform(&state[indx]);
      gitr_precision r2 = curand_uniform(&state[indx]);
      gitr_precision xsi = curand_uniform(&state[indx]);
#else
      std::normal_distribution<gitr_precision> distribution(0.0, 1.0);
      std::uniform_real_distribution<gitr_precision> dist(0.0, 1.0);
      gitr_precision n1 = distribution(state[indx]);
      gitr_precision n2 = distribution(state[indx]);
      gitr_precision r1 = dist(state[indx]);
      gitr_precision r2 = dist(state[indx]);
      gitr_precision r3 = dist(state[indx]);
      gitr_precision xsi = dist(state[indx]);
#endif

      gitr_precision ti_eV = cellTemperature[cellIndex];
      gitr_precision density = cellDensity[cellIndex];

      // get coulomb logarithm
      gitr_precision coulomb_log =
          CoulombLogarithm(density, ti_eV, particlesPointer->Z[indx], background_Z);

      // get relaxation time
      gitr_precision relaxation_time;
      relaxation_time =
          1.8e-19 * pow(background_amu * particlesPointer->amu[indx], 0.5) *
          pow(background_Z * particlesPointer->Z[indx], 2) * density * 1e-6 * coulomb_log /
          pow(ti_eV * background_amu + ti_eV * particlesPointer->amu[indx], 1.5);

      gitr_precision Prob = 1.0 - exp(-dt / relaxation_time);

      if (r1 <= Prob) {
        // calculate the center of mass velocity
        gitr_precision v_cm[3] = {0.0, 0.0, 0.0};
        v_cm[0] = (particlesPointer->amu[indx] * vx) / (particlesPointer->amu[indx] + background_amu);
        v_cm[1] = (particlesPointer->amu[indx] * vy) / (particlesPointer->amu[indx] + background_amu);
        v_cm[2] = (particlesPointer->amu[indx] * vz) / (particlesPointer->amu[indx] + background_amu);

        // get relative velocity
        gitr_precision v_rel[3] = {vx - v_cm[0], vy - v_cm[1], vz - v_cm[2]};
        gitr_precision v_rel_mag2 = v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2];
        gitr_precision v_rel_mag = sqrt(v_rel_mag2);

        // randomly select a particle in the same cell
        int randomIndx = getRandomParticleIndexInCell(cellIndex);

        // calculate the center of mass velocity of the selected particle
        gitr_precision v_cm_selected[3] = {0.0, 0.0, 0.0};
        v_cm_selected[0] = (particlesPointer->amu[randomIndx] * particlesPointer->vx[randomIndx]) /
                           (particlesPointer->amu[randomIndx] + background_amu);
        v_cm_selected[1] = (particlesPointer->amu[randomIndx] * particlesPointer->vy[randomIndx]) /
                           (particlesPointer->amu[randomIndx] + background_amu);
        v_cm_selected[2] = (particlesPointer->amu[randomIndx] * particlesPointer->vz[randomIndx]) /
                           (particlesPointer->amu[randomIndx] + background_amu);

        // get relative velocity with the selected particle
        gitr_precision v_rel_selected[3] = {particlesPointer->vx[randomIndx] - v_cm_selected[0],
                                            particlesPointer->vy[randomIndx] - v_cm_selected[1],
                                            particlesPointer->vz[randomIndx] - v_cm_selected[2]};
        gitr_precision v_rel_selected_mag2 =
            v_rel_selected[0] * v_rel_selected[0] + v_rel_selected[1] * v_rel_selected[1] +
            v_rel_selected[2] * v_rel_selected[2];
        gitr_precision v_rel_selected_mag = sqrt(v_rel_selected_mag2);

        // randomly scatter particle
        gitr_precision cos_theta = 2.0 * r2 - 1.0;
        gitr_precision sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        gitr_precision phi = 2.0 * M_PI * r3;
        gitr_precision v_scatter[3] = {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};

        // transform the scattering velocity to the center of mass frame
        for (int i = 0; i < 3; i++) {
          v_scatter[i] *= v_rel_selected_mag;
        }

        // transform the scattering velocity back to the lab frame and update the particle velocity
        particlesPointer->vx[indx] = vx + v_scatter[0];
        particlesPointer->vy[indx] = vy + v_scatter[1];
        particlesPointer->vz[indx] = vz + v_scatter[2];
      }
    }
  }
};
#endif
