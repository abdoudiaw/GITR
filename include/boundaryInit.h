#ifndef _BOUNDARYINIT_
#define _BOUNDARYINIT_

#include "pusher.h"
#include "interpolater.h"
#include "Boundary.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

// Constants

constexpr gitr_precision q_e = 1.60217662e-19;
constexpr gitr_precision me = 1.0/1836.152673440001;
constexpr gitr_precision PI_CONST = 3.141592653589793;
const gitr_precision VACUUM_PERMITTIVITY = 8.854187e-12;
const gitr_precision PI = PI_CONST;
const gitr_precision ELECTRON_MASS = me;
const gitr_precision CHARGE = q_e;

struct boundary_init {
    gitr_precision background_Z, background_amu;
    int nR_Temp, nZ_Temp, nx, nz, nxB, nzB, use_3d_geom;
    gitr_precision *ti, *te, *density, *ne, *bfieldR, *bfieldZ, *bfieldT;
    gitr_precision potential;


    // Simplified member initialization
    boundary_init(gitr_precision _background_Z, gitr_precision _background_amu,
                  gitr_precision* _density, gitr_precision* _ne, gitr_precision* _bfieldR,
                  gitr_precision* _bfieldZ, gitr_precision* _bfieldT,  
                  gitr_precision* _ti, gitr_precision* _te, 
                  gitr_precision _potential,
                  int _use_3d_geom)
        : background_Z{_background_Z}, background_amu{_background_amu}, ti{_ti}, te{_te}, potential{_potential}, use_3d_geom{_use_3d_geom},
          density{_density}, ne{_ne},  bfieldR{_bfieldR}, bfieldZ{_bfieldZ}, bfieldT{_bfieldT} {}

    void operator()(Boundary &b) const {
        gitr_precision midpointx;
        gitr_precision midpointy;
        gitr_precision midpointz;
    if( use_3d_geom )
        {
            midpointx = b.x1 + 0.666666667*(b.x2 + 0.5*(b.x3-b.x2)-b.x1);
            midpointy = b.y1 + 0.666666667*(b.y2 + 0.5*(b.y3-b.y2)-b.y1);
            midpointz = b.z1 + 0.666666667*(b.z2 + 0.5*(b.z3-b.z2)-b.z1);
        }
    else
        {
            midpointx = 0.5*(b.x2 - b.x1)+ b.x1;
            midpointy = 0.0;
            midpointz = 0.5*(b.z2 - b.z1) + b.z1;
        }
    b.density = density[0];
    b.ne = ne[0];
    b.ti =ti[0];
    b.te = te[0];
    gitr_precision B[3] = {0.0,0.0,0.0};

    B[0] = bfieldR[0];
    B[1] = bfieldT[0];
    B[2] = bfieldZ[0];

    gitr_precision magneticFieldNorm = vectorNorm(B);
    gitr_precision angleBetweenBAndSurfaceNormal = 0.0;

    if (use_3d_geom) {
        gitr_precision surfaceNormal[3] = {0.0, 0.0, 0.0};
        b.getSurfaceNormal(surfaceNormal, 0.0, 0.0, use_3d_geom);

        gitr_precision surfaceNormalNorm = vectorNorm(surfaceNormal);
        
        if (magneticFieldNorm == 0.0 || surfaceNormalNorm == 0.0) {
            angleBetweenBAndSurfaceNormal = 0.0;
            
            if (magneticFieldNorm == 0.0) {
                printf("Warning: Bfield is zero at surface %d\n", b.surfaceNumber);
            }
        } else {
            angleBetweenBAndSurfaceNormal = std::acos(vectorDotProduct(B, surfaceNormal) / (magneticFieldNorm * surfaceNormalNorm));
        }

        // Adjusting the angle if it's greater than pi/2
        if (angleBetweenBAndSurfaceNormal > PI / 2.0) {
            angleBetweenBAndSurfaceNormal = std::abs(angleBetweenBAndSurfaceNormal - PI);
        }

        b.unit_vec0 = b.inDir * b.a / b.plane_norm; 
        b.unit_vec1 = b.inDir * b.b / b.plane_norm; 
        b.unit_vec2 = b.inDir * b.c / b.plane_norm;
    }

    b.angle = angleBetweenBAndSurfaceNormal * 180.0 / PI;
    b.debyeLength = std::sqrt(VACUUM_PERMITTIVITY * b.te / (b.ne * std::pow(background_Z, 2) * CHARGE));

    gitr_precision sheathFactor = std::abs(0.5 * std::log((2 * PI * ELECTRON_MASS / background_amu) * (1 + b.ti / b.te)));
    gitr_precision normValue = std::acos(std::exp(-sheathFactor));

    if(b.ne == 0.0) b.debyeLength = 1.0e12;
    
    b.flux = 0.25*b.density*std::sqrt(8.0*b.ti*q_e/(PI_CONST*background_amu));

    b.impacts = 0.0;

    if(b.te == 0.0)
        b.potential = 0.0;
    else
        b.potential = sheathFactor*b.te;

    
    //  std::cout << "Surface number " << b.surfaceNumber << " has te "  << b.te << " and potential " << b.potential << " ne " << b.ne
    //     << " theta " << b.angle << std::endl;
    
    }	
};

#endif
