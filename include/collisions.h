#ifndef _COLLISIONS_
#define _COLLISIONS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "particles.h"
#include <cmath>
#include "interpolater.h"
#include "pusher.h"
#include "array.h"
#include "flags.h"
#include <cmath>
#include <limits>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#endif


gitr_precision CoulombLogarithm(gitr_precision ne, gitr_precision ti_eV, gitr_precision Z1, gitr_precision Z2){
  //note - this assumes ion and electron temperatures are the same
  
   if (ne == 0 || ti_eV == 0) {
    return 2;
  }

  double e2_cgs = 1.44e-7;       //eV-cm

  gitr_precision lam_e = 4 * M_PI * e2_cgs * ne * 1e-6/ ti_eV;            //1/cm^2, this is really 1/lam_e
  gitr_precision g_12 = Z1 * Z2 * e2_cgs * sqrt(lam_e) / ti_eV;

  gitr_precision K11;
  if (g_12 < 1.) {
    K11 = -0.25 * log(1.466 * g_12 - 1.7836 * pow(g_12, 2.0) + 1.4313 * pow(g_12, 3.0) - 0.55833 * pow(g_12, 4) + 0.061162 * pow(g_12, 5)) * 2.0;
  }
  else {
    K11 = (0.081033 - 0.091336 * log(g_12) + 0.05176 * pow(log(g_12), 2)) / (1.0 - 0.50026 * g_12 + 0.17044 * pow(g_12, 2)) * 2.0;
  }

  return fmax(K11, 2.0);
}


struct collisions {
    Particles* particlesPointer;
    gitr_precision dt;
    int nR_flowV, nY_flowV, nZ_flowV;
    gitr_precision* flowVGridr;
    gitr_precision* flowVGridy;
    gitr_precision* flowVGridz;
    gitr_precision* flowVr;
    gitr_precision* flowVz;
    gitr_precision* flowVt;
    int nR_Dens, nZ_Dens;
    gitr_precision* DensGridr;
    gitr_precision* DensGridz;
    gitr_precision* ni;
    int nR_Temp, nZ_Temp;
    gitr_precision* TempGridr;
    gitr_precision* TempGridz;
    gitr_precision* ti;
    gitr_precision* te;
    gitr_precision background_Z;
    gitr_precision background_amu;
    int nR_Bfield, nZ_Bfield;
    gitr_precision* BfieldGridR;
    gitr_precision* BfieldGridZ;
    gitr_precision* BfieldR;
    gitr_precision* BfieldZ;
    gitr_precision* BfieldT;
    Flags* gitr_flags;
    gitr_precision dv[3];

#ifdef __CUDACC__
    curandState* state;
#else
    std::mt19937* state;
#endif

    int flowv_interp;
    int field_aligned_values;

    collisions(Particles* _particlesPointer, gitr_precision _dt,
#ifdef __CUDACC__
                     curandState* _state,
#else
                     std::mt19937* _state,
#endif
                     int _nR_flowV, int _nY_flowV, int _nZ_flowV, gitr_precision* _flowVGridr,
                     gitr_precision* _flowVGridy, gitr_precision* _flowVGridz, gitr_precision* _flowVr,
                     gitr_precision* _flowVz, gitr_precision* _flowVt, int _nR_Dens, int _nZ_Dens,
                     gitr_precision* _DensGridr, gitr_precision* _DensGridz, gitr_precision* _ni,
                     int _nR_Temp, int _nZ_Temp, gitr_precision* _TempGridr, gitr_precision* _TempGridz,
                     gitr_precision* _ti, gitr_precision* _te, gitr_precision _background_Z,
                     gitr_precision _background_amu, int _nR_Bfield, int _nZ_Bfield,
                     gitr_precision* _BfieldGridR, gitr_precision* _BfieldGridZ,
                     gitr_precision* _BfieldR, gitr_precision* _BfieldZ, gitr_precision* _BfieldT,
                     Flags* _gitr_flags, int _flowv_interp,  int _field_aligned_values)
        : particlesPointer(_particlesPointer),
          dt(_dt),
          nR_flowV(_nR_flowV),
          nY_flowV(_nY_flowV),
          nZ_flowV(_nZ_flowV),
          flowVGridr(_flowVGridr),
          flowVGridy(_flowVGridy),
          flowVGridz(_flowVGridz),
          flowVr(_flowVr),
          flowVz(_flowVz),
          flowVt(_flowVt),
          nR_Dens(_nR_Dens),
          nZ_Dens(_nZ_Dens),
          DensGridr(_DensGridr),
          DensGridz(_DensGridz),
          ni(_ni),
          nR_Temp(_nR_Temp),
          nZ_Temp(_nZ_Temp),
          TempGridr(_TempGridr),
          TempGridz(_TempGridz),
          ti(_ti),
          te(_te),
          background_Z(_background_Z),
          background_amu(_background_amu),
          nR_Bfield(_nR_Bfield),
          nZ_Bfield(_nZ_Bfield),
          BfieldGridR(_BfieldGridR),
          BfieldGridZ(_BfieldGridZ),
          BfieldR(_BfieldR),
          BfieldZ(_BfieldZ),
          BfieldT(_BfieldT),
          gitr_flags(_gitr_flags),
          dv{0.0, 0.0, 0.0},
          state(_state),
          flowv_interp(_flowv_interp),
          field_aligned_values(_field_aligned_values) {}

    CUDA_CALLABLE_MEMBER_DEVICE
    void operator()(std::size_t indx) {
        if (particlesPointer->hitWall[indx] == 0.0 && particlesPointer->charge[indx] != 0.0) {

            gitr_precision flowVelocity[3] = {0.0};
            gitr_precision x = particlesPointer->xprevious[indx];
            gitr_precision y = particlesPointer->yprevious[indx];
            gitr_precision z = particlesPointer->zprevious[indx];
            gitr_precision vx = particlesPointer->vx[indx];
            gitr_precision vy = particlesPointer->vy[indx];
            gitr_precision vz = particlesPointer->vz[indx];

            if (flowv_interp == 3) {
                interp3dVector(&flowVelocity[0], x, y, z, nR_flowV, nY_flowV, nZ_flowV,
                               flowVGridr, flowVGridy, flowVGridz, flowVr, flowVz, flowVt);
            } else if (flowv_interp < 3) {
                if (field_aligned_values > 0) {
                    interpFieldAlignedVector(&flowVelocity[0], x, y, z, nR_flowV, nZ_flowV,
                                             flowVGridr, flowVGridz, flowVr,
                                             flowVz, flowVt, nR_Bfield, nZ_Bfield,
                                             BfieldGridR, BfieldGridZ, BfieldR,
                                             BfieldZ, BfieldT);
                } else {
                    interp2dVector(flowVelocity, x, y, z, nR_flowV, nZ_flowV,
                                   flowVGridr, flowVGridz, flowVr, flowVz, flowVt );
                }
            }

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
            gitr_precision ti_eV = interp2dCombined(x, y, z, nR_Temp, nZ_Temp, TempGridr, TempGridz, ti );
            gitr_precision density = interp2dCombined(x, y, z, nR_Dens, nZ_Dens, DensGridr, DensGridz, ni );

            // get coulomb logarithm
            gitr_precision coulomb_log = CoulombLogarithm( density, ti_eV, particlesPointer->Z[indx], background_Z);
     
            // get relaxation time
            gitr_precision relaxation_time;
            relaxation_time = 1.8e-19 * pow(background_amu * particlesPointer->amu[indx], 0.5) * pow(background_Z* particlesPointer->Z[indx], 2) * density * 1e-6 * coulomb_log / 
            pow(ti_eV * background_amu  + ti_eV * particlesPointer->amu[indx], 1.5) ;

            gitr_precision Prob = 1.0 - exp(-dt / relaxation_time);

            if (r1 <= Prob)
            {
                std::cout << "Collision occurred" << std::endl;
                // calculate the center of mass velocity
                gitr_precision v_cm[3] = {0.0, 0.0, 0.0};
                v_cm[0] = (particlesPointer->amu[indx] * vx + background_amu * flowVelocity[0]) / (particlesPointer->amu[indx] + background_amu);
                v_cm[1] = (particlesPointer->amu[indx] * vy + background_amu * flowVelocity[1]) / (particlesPointer->amu[indx] + background_amu);
                v_cm[2] = (particlesPointer->amu[indx] * vz + background_amu * flowVelocity[2]) / (particlesPointer->amu[indx] + background_amu);

                // get relative velocity
                gitr_precision v_rel= {0.0};
                gitr_precision v_rel_mag2 = (vx - v_cm[0]) * (vx - v_cm[0]) + (vy - v_cm[1]) * (vy - v_cm[1]) + (vz - v_cm[2]) * (vz - v_cm[2]);
                gitr_precision v_rel_mag = sqrt(v_rel_mag2);
    
                // randomly scatter particle
                gitr_precision cos_theta = 2.0 * r2 - 1.0;
                gitr_precision sin_theta = sqrt(1.0 - cos_theta * cos_theta);
                gitr_precision phi = 2.0 * M_PI * r3;
                gitr_precision vcm_scatter[3] = {  sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};

                // transform the scattering velocity to the center of mass frame
                gitr_precision v_scatter[3] = {0.0, 0.0, 0.0};
                for (int i = 0; i < 3; i++)
                {
                    v_scatter[i] *= v_rel_mag;
                }

                // transform the scattering velocity back to the lab frame and update the particle velocity
                particlesPointer->vx[indx] = vx + v_scatter[0];
                particlesPointer->vy[indx] = vy + v_scatter[1];
                particlesPointer->vz[indx] = vz + v_scatter[2];

                // print particle velocity
                std::cout << "vx = " << vx << std::endl;
                // print nR_Temp, nZ_Temp
                std::cout << "nR_Temp = " << nR_Temp << std::endl;
                std::cout << "nZ_Temp = " << nZ_Temp << std::endl;
            }
        }
    }
};
#endif
