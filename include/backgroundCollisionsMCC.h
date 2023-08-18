#ifndef _COULOMB_
#define _COULOMB_

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
#include <cmath>
#include <limits>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#endif
#include <fenv.h>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif



CUDA_CALLABLE_MEMBER
void getSlowDownFrequencies ( gitr_precision& nu_friction, gitr_precision& nu_deflection, gitr_precision& nu_parallel,
  gitr_precision& nu_energy, gitr_precision x, gitr_precision y,gitr_precision z, 
    gitr_precision vx, gitr_precision vy, gitr_precision vz, gitr_precision flowVr, gitr_precision flowVz,gitr_precision flowVt,
    gitr_precision charge, gitr_precision amu, gitr_precision ni, 
    gitr_precision ti, gitr_precision te, gitr_precision background_Z, gitr_precision background_amu,
    gitr_precision &T_background )
{
  //int feenableexcept(FE_INVALID | FE_OVERFLOW); //enables trapping of the floating-point exceptions
  gitr_precision Q = 1.60217662e-19;
  gitr_precision EPS0 = 8.854187e-12;
  gitr_precision MI = 1.6737236e-27;	
  gitr_precision ME = 9.10938356e-31;
        
  gitr_precision te_eV = te;
  gitr_precision ti_eV = ti;

  T_background = ti_eV;
  gitr_precision density = ni;

  gitr_precision flowVelocity[3]= {0.0};
  gitr_precision relativeVelocity[3] = {0.0, 0.0, 0.0};
  gitr_precision velocityNorm = 0.0;
  gitr_precision lam_d;
  gitr_precision lam;
  gitr_precision gam_electron_background;
  gitr_precision gam_ion_background;
  gitr_precision a_electron = 0.0;
  gitr_precision a_ion = 0.0;
  gitr_precision xx;
  gitr_precision psi_prime;
  gitr_precision psi_psiprime;
  gitr_precision psi;
  gitr_precision xx_e;
  gitr_precision psi_prime_e;
  gitr_precision psi_psiprime_e;
  gitr_precision psi_psiprime_psi2x = 0.0;
  gitr_precision psi_psiprime_psi2x_e = 0.0;
  gitr_precision psi_e;
  gitr_precision nu_0_i;
  gitr_precision nu_0_e;
  gitr_precision nu_friction_i;
  gitr_precision nu_deflection_i;
  gitr_precision nu_parallel_i;
  gitr_precision nu_energy_i;
  gitr_precision nu_friction_e;
  gitr_precision nu_deflection_e;
  gitr_precision nu_parallel_e;
  gitr_precision nu_energy_e;



                
  // charge  = 1.;
  // amu = 2.01410178;
  relativeVelocity[0] = vx - flowVr;
  relativeVelocity[1] = vy - flowVt;
  relativeVelocity[2] = vz - flowVz;

  velocityNorm = std::sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);                

  lam_d = std::sqrt(EPS0*te_eV/(density*std::pow(background_Z,2)*Q));//only one q in order to convert to J
  lam = 12.0*M_PI*density*std::pow(lam_d,3)/charge;
  gam_electron_background = 0.238762895*std::pow(charge,2)*std::log(lam)/(amu*amu);//constant = Q^4/(MI^2*4*pi*EPS0^2)
  gam_ion_background = 0.238762895*std::pow(charge,2)*std::pow(background_Z,2)*std::log(lam)/(amu*amu);//constant = Q^4/(MI^2*4*pi*EPS0^2)

  if(gam_electron_background < 0.0) gam_electron_background=0.0;
  if(gam_ion_background < 0.0) gam_ion_background=0.0;
  a_ion = background_amu*MI/(2*ti_eV*Q);// %q is just to convert units - no z needed
  a_electron = ME/(2*te_eV*Q);// %q is just to convert units - no z needed

  xx = std::pow(velocityNorm,2)*a_ion;
  psi_prime = 2.0*std::sqrt(xx/M_PI)*std::exp(-xx);
  psi_psiprime = std::erf(std::sqrt(xx));
  psi = psi_psiprime - psi_prime;
  
  xx_e = std::pow(velocityNorm,2)*a_electron;
  psi_prime_e = 1.128379*std::sqrt(xx_e);
  psi_e = 0.75225278*std::pow(xx_e,1.5);
  psi_psiprime_e = psi_e+psi_prime_e;
  psi_psiprime_psi2x_e = 1.128379*std::sqrt(xx_e)*expf(-xx_e);
  
  nu_0_i = gam_electron_background*density/std::pow(velocityNorm,3);
  nu_0_e = gam_ion_background*density/std::pow(velocityNorm,3);
  
  nu_friction_i = (1+amu/background_amu)*psi*nu_0_i;
  nu_deflection_i = 2*(psi_psiprime - psi/(2*xx))*nu_0_i;
  nu_parallel_i = psi/xx*nu_0_i;
  nu_energy_i = 2*(amu/background_amu*psi - psi_prime)*nu_0_i;

  
  nu_friction_e = (1+amu/(ME/MI))*psi_e*nu_0_e;
  nu_deflection_e = 2*(psi_psiprime_psi2x_e)*nu_0_e;
  nu_parallel_e = psi_e/xx_e*nu_0_e;
  nu_energy_e = 2*(amu/(ME/MI)*psi_e - psi_prime_e)*nu_0_e;
                    
  nu_friction = nu_friction_i ;//+ nu_friction_e;
  nu_deflection = nu_deflection_i ;//+ nu_deflection_e;
  nu_parallel = nu_parallel_i;// + nu_parallel_e;
  nu_energy = nu_energy_i;// + nu_energy_e;
    
  if(te_eV <= 0.0 || ti_eV <= 0.0)
  {
    nu_friction = 0.0;
    nu_deflection = 0.0;
    nu_parallel = 0.0;
    nu_energy = 0.0;
  }
  if(density <= 0.0)
  {
    nu_friction = 0.0;
    nu_deflection = 0.0;
    nu_parallel = 0.0;
    nu_energy = 0.0;
  }
}
CUDA_CALLABLE_MEMBER
void getSlowDownDirections2 (gitr_precision parallel_direction[], gitr_precision perp_direction1[], gitr_precision perp_direction2[],
        gitr_precision vx, gitr_precision vy, gitr_precision vz)
{
  gitr_precision v = std::sqrt(vx*vx + vy*vy + vz*vz);
  
  if(v == 0.0)
  {
    v = 1.0;
    vz = 1.0;
    vx = 0.0;
    vy = 0.0;
  }
  
  gitr_precision ez1 = vx/v;
  gitr_precision ez2 = vy/v;
  gitr_precision ez3 = vz/v;
    
  // Get perpendicular velocity unit vectors
  // this comes from a cross product of
  // (ez1,ez2,ez3)x(0,0,1)
  gitr_precision ex1 = ez2;
  gitr_precision ex2 = -ez1;
  gitr_precision ex3 = 0.0;
  
  // The above cross product will be zero for particles
  // with a pure z-directed (ez3) velocity
  // here we find those particles and get the perpendicular 
  // unit vectors by taking the cross product
  // (ez1,ez2,ez3)x(0,1,0) instead
  gitr_precision exnorm = std::sqrt(ex1*ex1 + ex2*ex2);
  if(std::abs(exnorm) < 1.0e-12){
  ex1 = -ez3;
  ex2 = 0.0;
  ex3 = ez1;
  }
  // Ensure all the perpendicular direction vectors
  // ex are unit
  exnorm = std::sqrt(ex1*ex1+ex2*ex2 + ex3*ex3);
  ex1 = ex1/exnorm;
  ex2 = ex2/exnorm;
  ex3 = ex3/exnorm;
  
  // Find the second perpendicular direction 
  // by taking the cross product
  // (ez1,ez2,ez3)x(ex1,ex2,ex3)
  gitr_precision ey1 = ez2*ex3 - ez3*ex2;
  gitr_precision ey2 = ez3*ex1 - ez1*ex3;
  gitr_precision ey3 = ez1*ex2 - ez2*ex1;
  parallel_direction[0] = ez1; 
  parallel_direction[1] = ez2;
  parallel_direction[2] = ez3;
  
  perp_direction1[0] = ex1; 
  perp_direction1[1] = ex2;
  perp_direction1[2] = ex3;
  
  perp_direction2[0] = ey1; 
  perp_direction2[1] = ey2;
  perp_direction2[2] = ey3;
}


struct backgroundCollisionsMCC { 
    Particles *particlesPointer;
    gitr_precision flowVr;
    gitr_precision flowVz;
    gitr_precision flowVt;
    gitr_precision ne;
    gitr_precision ni;
    gitr_precision ti;
    gitr_precision te;
    gitr_precision background_Z;
    gitr_precision background_amu;
    gitr_precision dv[3];
    #if __CUDACC__
            curandState *state;
    #else
            std::mt19937 *state;
    #endif

    backgroundCollisionsMCC(Particles *_particlesPointer,
    #if __CUDACC__
                                curandState *_state,
    #else
                                std::mt19937 *_state,
    #endif
            gitr_precision _flowVr, gitr_precision _flowVz, gitr_precision _flowVt,
            gitr_precision _ne,
            gitr_precision _ni,
            gitr_precision _ti, gitr_precision _te,
            gitr_precision _background_Z, gitr_precision _background_amu )
          : particlesPointer(_particlesPointer),
            flowVr(_flowVr),
            flowVz(_flowVz),
            flowVt(_flowVt),
            ne(_ne),
            ni(_ni),
            ti(_ti),
            te(_te),
            background_Z(_background_Z),
            background_amu(_background_amu),
            dv{0.0, 0.0, 0.0},
            state(_state)
            { }

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) { 

  if(particlesPointer->hitWall[indx] == 0.0 && particlesPointer->charge[indx] != 0.0)
  {

    gitr_precision dt = particlesPointer->dt[indx];	   
    gitr_precision T_background = 0.0;
    gitr_precision nu_friction = 0.0;
    gitr_precision nu_deflection = 0.0;
    gitr_precision nu_parallel = 0.0;
    gitr_precision nu_energy = 0.0;
    gitr_precision flowVelocity[3]= {0.0};
    gitr_precision relativeVelocity[3] = {0.0};
    gitr_precision velocityCollisions[3]= {0.0};	
    gitr_precision velocityRelativeNorm;	
    gitr_precision parallel_direction[3] = {0.0};
    gitr_precision perp_direction1[3] = {0.0};
    gitr_precision perp_direction2[3] = {0.0};

    
    gitr_precision x = particlesPointer->xprevious[indx];
    gitr_precision y = particlesPointer->yprevious[indx];
    gitr_precision z = particlesPointer->zprevious[indx];
    gitr_precision vx = particlesPointer->vx[indx];
    gitr_precision vy = particlesPointer->vy[indx];
    gitr_precision vz = particlesPointer->vz[indx];

    relativeVelocity[0] = vx - flowVr;
    relativeVelocity[1] = vy - flowVt;
    relativeVelocity[2] = vz - flowVz;

    velocityRelativeNorm = vectorNorm(relativeVelocity);

#ifdef __CUDACC__
    gitr_precision n1 = curand_normal(&state[indx]);
    gitr_precision n2 = curand_normal(&state[indx]);
    gitr_precision r1 = curand_uniform(&state[indx]);
    gitr_precision r2 = curand_uniform(&state[indx]);
    gitr_precision r3 = curand_uniform(&state[indx]);
    gitr_precision xsi = curand_uniform(&state[indx]);
#else
    std::normal_distribution<gitr_precision> distribution(0.0,1.0);
    std::uniform_real_distribution<gitr_precision> dist(0.0, 1.0);
    gitr_precision n1 = distribution(state[indx]);
    gitr_precision n2 = distribution(state[indx]);
    gitr_precision r1 = dist(state[indx]);
    gitr_precision r2 = dist(state[indx]);
    gitr_precision r3 = dist(state[indx]);
    gitr_precision xsi = dist(state[indx]);
#endif

    getSlowDownFrequencies(nu_friction, nu_deflection, nu_parallel, nu_energy,
                             x, y, z, vx, vy, vz,
                            flowVr, flowVz, flowVt, particlesPointer->charge[indx], particlesPointer->amu[indx], ni, ti, te, background_Z, background_amu,
                             T_background );

    getSlowDownDirections2(parallel_direction, perp_direction1, perp_direction2,
                            relativeVelocity[0] , relativeVelocity[1] , relativeVelocity[2] );
    
{

    gitr_precision ti_eV = ti;
    gitr_precision density = ni;

    
    if(nu_parallel <=0.0) nu_parallel = 0.0;

    gitr_precision coeff_par = n1 * std::sqrt(2.0*nu_parallel * dt);
    gitr_precision cosXsi = cos(2.0 * M_PI * xsi) - 0.0028;

    if(cosXsi > 1.0) cosXsi = 1.0;
    gitr_precision sinXsi = sin(2.0 * M_PI * xsi);

    if(nu_deflection <=0.0) nu_deflection = 0.0;

    gitr_precision coeff_perp1 = cosXsi * std::sqrt(nu_deflection * dt*0.5);
    gitr_precision coeff_perp2 = sinXsi * std::sqrt(nu_deflection * dt*0.5);
      
    gitr_precision nuEdt = nu_energy * dt;
    if (nuEdt < -1.0) nuEdt = -1.0;

    gitr_precision vx_relative = velocityRelativeNorm*(1.0-0.5*nuEdt)*((1.0 + coeff_par) * parallel_direction[0] + std::abs(n2)*(coeff_perp1 * perp_direction1[0] + coeff_perp2 * perp_direction2[0])) - velocityRelativeNorm*dt*nu_friction*parallel_direction[0];
    gitr_precision vy_relative = velocityRelativeNorm*(1.0-0.5*nuEdt)*((1.0 + coeff_par) * parallel_direction[1] + std::abs(n2)*(coeff_perp1 * perp_direction1[1] + coeff_perp2 * perp_direction2[1])) - velocityRelativeNorm*dt*nu_friction*parallel_direction[1];
    gitr_precision vz_relative = velocityRelativeNorm*(1.0-0.5*nuEdt)*((1.0 + coeff_par) * parallel_direction[2] + std::abs(n2)*(coeff_perp1 * perp_direction1[2] + coeff_perp2 * perp_direction2[2])) - velocityRelativeNorm*dt*nu_friction*parallel_direction[2];

    particlesPointer->vx[indx] = vx_relative + flowVelocity[0]; 
    particlesPointer->vy[indx] = vy_relative + flowVelocity[1]; 
    particlesPointer->vz[indx] = vz_relative + flowVelocity[2];

    this->dv[0] = velocityCollisions[0];
    this->dv[1] = velocityCollisions[1];
    this->dv[2] = velocityCollisions[2];

    particlesPointer->x[indx] = x +  particlesPointer->vx[indx]*dt;
    particlesPointer->y[indx] = y +  particlesPointer->vy[indx]*dt;
    particlesPointer->z[indx] = z +  particlesPointer->vz[indx]*dt;

    }
  }
}
};

#endif
