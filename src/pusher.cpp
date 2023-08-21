/* For better separation and code readibility,
TODO: have a separate function to calculate total electric field
 Same for the magnetic field
 Aggregate total forces into F_accelerating ~E and gyro F_~B
 Pass that to pusher/Boris for particles update. 
 */

#include <vector>
#include <algorithm>
#include <cassert>

#include <tuple>
#include <cmath>

#include "pusher.h"
#include "sheathModel.h"
#include "constants.h"


// Compute cross product of two vectors: AxB
CUDA_CALLABLE_MEMBER
gitr_precision compute_cross_product(const gitr_precision A[3], const gitr_precision B[3], int i) {
    return A[(i + 1) % 3] * B[(i + 2) % 3] - A[(i + 2) % 3] * B[(i + 1) % 3];
}

CUDA_CALLABLE_MEMBER
void vectorAdd(gitr_precision A[], gitr_precision B[],gitr_precision C[])
{
    C[0] = A[0] + B[0];
    C[1] = A[1] + B[1];
    C[2] = A[2] + B[2];
}

CUDA_CALLABLE_MEMBER
void vectorSubtract(gitr_precision A[], gitr_precision B[],gitr_precision C[])
{
    C[0] = A[0] - B[0];
    C[1] = A[1] - B[1];
    C[2] = A[2] - B[2];
}

CUDA_CALLABLE_MEMBER
void vectorScalarMult(gitr_precision a, gitr_precision B[],gitr_precision C[])
{
    C[0] = a*B[0];
    C[1] = a*B[1];
    C[2] = a*B[2];
}

CUDA_CALLABLE_MEMBER
void vectorAssign(gitr_precision a, gitr_precision b,gitr_precision c, gitr_precision D[])
{
    D[0] = a;
    D[1] = b;
    D[2] = c;
}

CUDA_CALLABLE_MEMBER
gitr_precision vectorNorm(gitr_precision A[])
{
    gitr_precision norm = 0.0;
    norm = std::sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);

    return norm;
}

CUDA_CALLABLE_MEMBER
void vectorNormalize(gitr_precision A[],gitr_precision B[])
{
    gitr_precision norm = 0.0;
    norm = std::sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    B[0] = A[0]/norm;
    B[1] = A[1]/norm;
    B[2] = A[2]/norm;

}

CUDA_CALLABLE_MEMBER
gitr_precision vectorDotProduct(gitr_precision A[], gitr_precision B[])
{
    gitr_precision c = A[0]*B[0] +  A[1]*B[1] + A[2]*B[2];
    return c;
}

CUDA_CALLABLE_MEMBER
void vectorCrossProduct(gitr_precision A[], gitr_precision B[], gitr_precision C[])
{
    gitr_precision tmp[3] = {0.0,0.0,0.0};
    tmp[0] = A[1]*B[2] - A[2]*B[1];
    tmp[1] = A[2]*B[0] - A[0]*B[2];
    tmp[2] = A[0]*B[1] - A[1]*B[0];

    C[0] = tmp[0];
    C[1] = tmp[1];
    C[2] = tmp[2];
}


CUDA_CALLABLE_MEMBER
gitr_precision getSheathElectricField ( gitr_precision x0, gitr_precision y, gitr_precision z, gitr_precision E[], Boundary *boundaryVector, int nLines,
         int&  closestBoundaryIndex, int use_3d_geom, gitr_precision mass, gitr_precision charge, gitr_precision Bmagnitude)
    {

    gitr_precision pot = 0.0;
    int minIndex = 0;
    gitr_precision minDistance = 1.0e12;
    gitr_precision Emag = 0.0;
    gitr_precision angle = 0.0;
    gitr_precision fd = 0.0;
    gitr_precision directionUnitVector[ 3 ] = { 0.0, 0.0, 0.0 };
    gitr_precision Er = 0.0;
    gitr_precision Et = 0.0;

    if( use_3d_geom > 0 )
    {
    gitr_precision p0[3] = {x0,y,z};
      gitr_precision a = 0.0;
      gitr_precision b = 0.0;
      gitr_precision c = 0.0;
      gitr_precision d = 0.0;
      gitr_precision plane_norm = 0.0;
      gitr_precision pointToPlaneDistance0 = 0.0;
      gitr_precision pointToPlaneDistance1 = 0.0;
      gitr_precision signPoint0 = 0.0;
      gitr_precision signPoint1 = 0.0;
      gitr_precision t = 0.0;
      gitr_precision A[3] = {0.0,0.0,0.0};
      gitr_precision B[3] = {0.0,0.0,0.0};
      gitr_precision C[3] = {0.0,0.0,0.0};
      gitr_precision AB[3] = {0.0,0.0,0.0};
      gitr_precision AC[3] = {0.0,0.0,0.0};
      gitr_precision BC[3] = {0.0,0.0,0.0};
      gitr_precision CA[3] = {0.0,0.0,0.0};
      gitr_precision p[3] = {0.0,0.0,0.0};
      gitr_precision Ap[3] = {0.0,0.0,0.0};
      gitr_precision Bp[3] = {0.0,0.0,0.0};
      gitr_precision Cp[3] = {0.0,0.0,0.0};
      gitr_precision p0A[3] = {0.0,0.0,0.0};
      gitr_precision p0B[3] = {0.0,0.0,0.0};
      gitr_precision p0C[3] = {0.0,0.0,0.0};
      gitr_precision p0AB[3] = {0.0,0.0,0.0};
      gitr_precision p0BC[3] = {0.0,0.0,0.0};
      gitr_precision p0CA[3] = {0.0,0.0,0.0};
      gitr_precision p0Anorm = 0.0;
      gitr_precision p0Bnorm = 0.0;
      gitr_precision p0Cnorm = 0.0;
      gitr_precision normalVector[3] = {0.0,0.0,0.0};
      gitr_precision crossABAp[3] = {0.0,0.0,0.0};
      gitr_precision crossBCBp[3] = {0.0,0.0,0.0};
      gitr_precision crossCACp[3] = {0.0,0.0,0.0};
      gitr_precision dot0 = 0.0;
      gitr_precision dot1 = 0.0;
      gitr_precision dot2 = 0.0;

      gitr_precision normAB = 0.0;
      gitr_precision normBC = 0.0;
      gitr_precision normCA = 0.0;
      gitr_precision ABhat[3] = {0.0,0.0,0.0};
      gitr_precision BChat[3] = {0.0,0.0,0.0};
      gitr_precision CAhat[3] = {0.0,0.0,0.0};
      gitr_precision tAB = 0.0;
      gitr_precision tBC = 0.0;
      gitr_precision tCA = 0.0;
      gitr_precision projP0AB[3] = {0.0,0.0,0.0};
      gitr_precision projP0BC[3] = {0.0,0.0,0.0};
      gitr_precision projP0CA[3] = {0.0,0.0,0.0};
      gitr_precision p0ABdist = 0.0;
      gitr_precision p0BCdist = 0.0;
      gitr_precision p0CAdist = 0.0;
      gitr_precision perpDist = 0.0;
      gitr_precision signDot0 = 0.0;
      gitr_precision signDot1 = 0.0;
      gitr_precision signDot2 = 0.0;
      gitr_precision totalSigns = 0.0;
      minDistance = 1.0e12;
      int nBoundariesCrossed = 0;
      int boundariesCrossed[6] = {0,0,0,0,0,0};
      gitr_precision distances[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      gitr_precision normals[21] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                           0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                           0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  int top_limit = -1;

  gitr_precision dr;
  gitr_precision dy;
  gitr_precision dz;

  int rInd;
  int yInd;
  int zInd;
top_limit = nLines;

  for (int k=0; k < top_limit; k++) //n_closeGeomElements
  {
    int i = -1;

      i = k;
    

    a = boundaryVector[i].a;
    b = boundaryVector[i].b;
    c = boundaryVector[i].c;
    d = boundaryVector[i].d;
    plane_norm = boundaryVector[i].plane_norm;
    pointToPlaneDistance0 = (a * p0[0] + b * p0[1] + c * p0[2] + d) / plane_norm;
    vectorAssign(a / plane_norm, b / plane_norm, c / plane_norm, normalVector);
    vectorAssign(p0[0] - pointToPlaneDistance0 * normalVector[0],
                 p0[1] - pointToPlaneDistance0 * normalVector[1],
                 p0[2] - pointToPlaneDistance0 * normalVector[2], p);

    vectorAssign(boundaryVector[i].x1, boundaryVector[i].y1,
                 boundaryVector[i].z1, A);
    vectorAssign(boundaryVector[i].x2, boundaryVector[i].y2,
                 boundaryVector[i].z2, B);
    vectorAssign(boundaryVector[i].x3, boundaryVector[i].y3,
                 boundaryVector[i].z3, C);

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
    signDot0 = std::copysign(1.0,vectorDotProduct(crossABAp, normalVector));
    signDot1 = std::copysign(1.0,vectorDotProduct(crossBCBp, normalVector));
    signDot2 = std::copysign(1.0,vectorDotProduct(crossCACp, normalVector));
         totalSigns = std::abs(signDot0 + signDot1 + signDot2);
         vectorSubtract(A,p0,p0A);
         vectorSubtract(B,p0,p0B);
         vectorSubtract(C,p0,p0C);
         
         p0Anorm = vectorNorm(p0A);   
         p0Bnorm = vectorNorm(p0B);   
         p0Cnorm = vectorNorm(p0C);
         distances[1] = p0Anorm;   
         distances[2] = p0Bnorm;   
         distances[3] = p0Cnorm;   
             normals[3] = p0A[0]/p0Anorm;
             normals[4] = p0A[1]/p0Anorm;
             normals[5] = p0A[2]/p0Anorm;
             normals[6] = p0B[0]/p0Bnorm;
             normals[7] = p0B[1]/p0Bnorm;
             normals[8] = p0B[2]/p0Bnorm;
             normals[9] = p0C[0]/p0Cnorm;
             normals[10] = p0C[1]/p0Cnorm;
             normals[11] = p0C[2]/p0Cnorm;
         normAB = vectorNorm(AB);
         normBC = vectorNorm(BC);
         normCA = vectorNorm(CA);
         vectorAssign(AB[0]/normAB,AB[1]/normAB,AB[2]/normAB,ABhat);
         vectorAssign(BC[0]/normBC,BC[1]/normBC,BC[2]/normBC,BChat);
         vectorAssign(CA[0]/normCA,CA[1]/normCA,CA[2]/normCA,CAhat);
         
         tAB = vectorDotProduct(p0A,ABhat);
         tBC = vectorDotProduct(p0B,BChat);
         tCA = vectorDotProduct(p0C,CAhat);
         tAB = -1.0*tAB;
         tBC = -1.0*tBC;
         tCA = -1.0*tCA;
         if((tAB > 0.0) && (tAB < normAB))
         {
             vectorScalarMult(tAB,ABhat,projP0AB);
             vectorAdd(A,projP0AB,projP0AB);
             vectorSubtract(projP0AB,p0,p0AB);
             p0ABdist = vectorNorm(p0AB);
             distances[4] = p0ABdist;   
             normals[12] = p0AB[0]/p0ABdist;
             normals[13] = p0AB[1]/p0ABdist;
             normals[14] = p0AB[2]/p0ABdist;

         }
         else
         {
             p0ABdist = 1.0e12;
             distances[4] = p0ABdist;   
         } 
         
         
         if((tBC > 0.0) && (tBC < normBC))
         {
             vectorScalarMult(tBC,ABhat,projP0BC);
             vectorAdd(B,projP0BC,projP0BC);
             vectorSubtract(projP0BC,p0,p0BC);
             p0BCdist = vectorNorm(p0BC);
             distances[5] = p0BCdist;   
             normals[15] = p0BC[0]/p0BCdist;
             normals[16] = p0BC[1]/p0BCdist;
             normals[17] = p0BC[2]/p0BCdist;

         }
         else
         {
             p0BCdist = 1.0e12;
             distances[5] = p0BCdist;   

         } 
         
         if((tCA > 0.0) && (tCA < normCA))
         {
             vectorScalarMult(tCA,CAhat,projP0CA);
             vectorAdd(C,projP0CA,projP0CA);
             vectorSubtract(projP0CA,p0,p0CA);
             p0CAdist = vectorNorm(p0CA);
             distances[6] = p0CAdist;   
             normals[18] = p0CA[0]/p0CAdist;
             normals[19] = p0CA[1]/p0CAdist;
             normals[20] = p0CA[2]/p0CAdist;
         }
         else
         {
             p0CAdist = 1.0e12;
             distances[6] = p0CAdist;   
         } 

         if (totalSigns == 3.0)
         {
                perpDist = std::abs(pointToPlaneDistance0); 
                vectorSubtract(p,p0 ,normalVector);
                vectorNormalize(normalVector,normalVector);
             distances[0] = perpDist;   
             normals[0] = boundaryVector[i].unit_vec0; //normalVector[0];
             normals[1] = boundaryVector[i].unit_vec1; //normalVector[1];
             normals[2] = boundaryVector[i].unit_vec2; //normalVector[2];
             //}
         }
         else
         {
             perpDist = 1.0e12;
             distances[0] = perpDist;   
         }
         int index = 0;
         for(int j = 0; j < 7; j++)
         {
            if(distances[j] < distances[index])
            index = j;              
         }

         if (distances[index] < minDistance)
         {
                 minDistance = distances[index];
                 vectorAssign(normals[index*3], normals[index*3+1],normals[index*3+2], directionUnitVector);
          closestBoundaryIndex = i;
          minIndex = i;
         }
  }

    if(isnan(directionUnitVector[0]) || isnan(directionUnitVector[1]) || isnan(directionUnitVector[2])){
	    directionUnitVector[0] = 0.0;
	    directionUnitVector[1] = 0.0;
	    directionUnitVector[2] = 0.0;
    }
    }
    else
    {
                
    int direction_type;
    gitr_precision tol = 1e12;
    gitr_precision point1_dist;
    gitr_precision point2_dist;
    gitr_precision perp_dist;
    gitr_precision vectorMagnitude;
    gitr_precision max = 0.0;
    gitr_precision min = 0.0;
    gitr_precision Bfabsfperp = 0.0;
    gitr_precision distanceToParticle = 0.0;
    int pointLine=0;
    gitr_precision x;
    x = x0;
    
    int top_limit = -1;
    gitr_precision dr;
    gitr_precision dz;

    int rInd;
    int zInd;

  top_limit = nLines;
  
  for( int k = 0; k < top_limit; k++) //n_closeGeomElements
    {
      int j = -1;

     j = k;

       gitr_precision boundZhere = boundaryVector[j].Z;
       
        if (boundZhere != 0.0)
        {
            point1_dist = std::sqrt((x - boundaryVector[j].x1)*(x - boundaryVector[j].x1) + 
                    (z - boundaryVector[j].z1)*(z - boundaryVector[j].z1));
            point2_dist = std::sqrt((x - boundaryVector[j].x2)*(x - boundaryVector[j].x2) + 
                                        (z - boundaryVector[j].z2)*(z - boundaryVector[j].z2));
            perp_dist = (boundaryVector[j].slope_dzdx*x - z + boundaryVector[j].intercept_z)/
                std::sqrt(boundaryVector[j].slope_dzdx*boundaryVector[j].slope_dzdx + 1.0);   
	
	
          if (std::abs(boundaryVector[j].slope_dzdx) >= tol*0.75)
	  {
	   perp_dist = x0 - boundaryVector[j].x1;
	  }
            if (point1_dist > point2_dist)
            {
                max = point1_dist;
                min = point2_dist;
            }
            else
            {
                max = point2_dist;
                min = point1_dist;
            }
            if (boundaryVector[j].length*boundaryVector[j].length + perp_dist*perp_dist >=
                    max*max)
            {
                distanceToParticle = std::abs(perp_dist);
                pointLine = 1;
            }
            else
            {
                distanceToParticle = min;
                if (boundaryVector[j].distanceToParticle == point1_dist)
                {
                    pointLine = 2;
                }
                else
                {
                    pointLine = 3;
                }
            }

            if (distanceToParticle < minDistance)
            {
                minDistance = distanceToParticle;
                minIndex = j;
                closestBoundaryIndex = j;
                direction_type = pointLine;
            }
        }
        else
        {
            distanceToParticle = tol;
        }
    }
    if (direction_type == 1)
    {
        if (boundaryVector[minIndex].slope_dzdx == 0)
        {
            directionUnitVector[0] = 0.0;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 1.0 * std::copysign(1.0,boundaryVector[minIndex].z1 - z);
        }
        else if (std::abs(boundaryVector[minIndex].slope_dzdx)>= 0.75*tol)
        {
            
            directionUnitVector[0] = boundaryVector[minIndex].x1 - x;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 0.0;
        }
        else
        {
            directionUnitVector[0] = 1.0 * std::copysign(1.0,(z - boundaryVector[minIndex].intercept_z)/(boundaryVector[minIndex].slope_dzdx) - x0);
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 1.0 * std::copysign(1.0,perp_dist)/(boundaryVector[minIndex].slope_dzdx);
        }
    }
    else if (direction_type == 2)
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x1 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z1 - z);
    }
    else
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x2 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z2 - z);
    }

    vectorMagnitude = std::sqrt(directionUnitVector[0]*directionUnitVector[0] + directionUnitVector[1]*directionUnitVector[1]
                                + directionUnitVector[2]*directionUnitVector[2]);
    directionUnitVector[0] = directionUnitVector[0]/vectorMagnitude;
    directionUnitVector[1] = directionUnitVector[1]/vectorMagnitude;
    directionUnitVector[2] = directionUnitVector[2]/vectorMagnitude;
    }
    if(minDistance == 0.0 || boundaryVector[minIndex].ne == 0.0 || boundaryVector[minIndex].te == 0.0)
        {
            Emag = 0.0;
            directionUnitVector[0] = 0.0;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 0.0;
        }
    else {
        // Call sheath model
        gitr_precision phiNormalizationFactor =  11600.0*gitr_constants::boltz/boundaryVector[minIndex].te/gitr_constants::p_mass;

        gitr_precision omega_c = charge * gitr_constants::e_charge * Bmagnitude / mass / gitr_constants::p_mass;
        gitr_precision cs = sqrt(gitr_constants::boltz * boundaryVector[minIndex].te * 11600. / mass / gitr_constants::p_mass);
        gitr_precision larmorRadius =  cs / omega_c;
        gitr_precision E_new = BrooksModel(abs(minDistance), boundaryVector[minIndex].angle, larmorRadius, boundaryVector[minIndex].te);
        Emag = abs(E_new);
        // printf ("Emag: %g \n", Emag);

    }
    Er = Emag*directionUnitVector[0];
    Et = Emag*directionUnitVector[1];
    E[2] = Emag*directionUnitVector[2];

    if( use_3d_geom > 0 )
    {
            E[0] = Er;
            E[1] = Et;

    }
    else
    {

            E[0] = Er;
            E[1] = Et;

    }

      return minDistance;
}

pusher::pusher(
  Particles *_particlesPointer,
  Boundary *_boundaryVector,
  int _nLines,
  int _nR_Bfield,
  int _nZ_Bfield,
  gitr_precision * _BfieldGridRDevicePointer,
  gitr_precision * _BfieldGridZDevicePointer,
  gitr_precision * _BfieldRDevicePointer,
  gitr_precision * _BfieldZDevicePointer,
  gitr_precision * _BfieldTDevicePointer,
 Flags* _gitr_flags)
  : 
  particlesPointer(_particlesPointer),
        boundaryVector(_boundaryVector),
        nLines(_nLines),
        nR_Bfield(_nR_Bfield),
        nZ_Bfield(_nZ_Bfield),
        BfieldGridRDevicePointer(_BfieldGridRDevicePointer),
        BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
        BfieldRDevicePointer(_BfieldRDevicePointer),
        BfieldZDevicePointer(_BfieldZDevicePointer),
        BfieldTDevicePointer(_BfieldTDevicePointer),
       gitr_flags(_gitr_flags)
        {}

CUDA_CALLABLE_MEMBER    
void pusher::operator()(std::size_t indx) {
    gitr_precision position[3] = {particlesPointer->xprevious[indx], particlesPointer->yprevious[indx], particlesPointer->zprevious[indx]};

    gitr_precision v[3] = {particlesPointer->vx[indx], particlesPointer->vy[indx], particlesPointer->vz[indx]};
    gitr_precision E[3] = {0.0, 0.0, 0.0};
    gitr_precision B[3] = {0.0, 0.0, 0.0};

    int closestBoundaryIndex;
    int geomIndex = 0;
  
    // if neutral stream 
    if (particlesPointer->charge[indx] == 0.0) {
        particlesPointer->x[indx] += particlesPointer->vx[indx] * particlesPointer->dt[indx];
        particlesPointer->y[indx] += particlesPointer->vy[indx] * particlesPointer->dt[indx];
        particlesPointer->z[indx] += particlesPointer->vz[indx] * particlesPointer->dt[indx];
        return;
    }
     // Get magnetic field
    interp2dVector(&B[0], position[0], position[1], position[2], nR_Bfield, nZ_Bfield, BfieldGridRDevicePointer, BfieldGridZDevicePointer, BfieldRDevicePointer, BfieldZDevicePointer, BfieldTDevicePointer );
    gitr_precision q_over_m = particlesPointer->charge[indx] * gitr_constants::e_charge / particlesPointer->amu[indx] / gitr_constants::p_mass;
    gitr_precision dt = particlesPointer->dt[indx];

   // get larmorRadius
    gitr_precision Bmagnitude = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

    // get electric field
    if (gitr_flags->USE_SHEATHEFIELD ) {
        getSheathElectricField(position[0], position[1], position[2], E, boundaryVector, nLines, closestBoundaryIndex, gitr_flags->USE3DTETGEOM, particlesPointer->amu[indx], particlesPointer->charge[indx], Bmagnitude );
    }

    // call the boris_push 
    if (particlesPointer->hitWall[indx] == 0.0)
    {
        std::tie(particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx], particlesPointer->vx[indx], particlesPointer->vy[indx], particlesPointer->vz[indx]) = 
        borisPush(particlesPointer, indx, E, B);
    }
}

std::tuple<gitr_precision, gitr_precision, gitr_precision,
           gitr_precision, gitr_precision, gitr_precision> 
borisPush(Particles *particlesPointer, std::size_t indx,
          gitr_precision E[3], gitr_precision B[3]) {

    // Extract initial values from particle using pointers and indx
    gitr_precision x_i = particlesPointer->x[indx];
    gitr_precision y_i = particlesPointer->y[indx];
    gitr_precision z_i = particlesPointer->z[indx];
    gitr_precision ux_i = particlesPointer->vx[indx];
    gitr_precision uy_i = particlesPointer->vy[indx];
    gitr_precision uz_i = particlesPointer->vz[indx];

    // Extract magnetic field components
    gitr_precision Bx = B[0];
    gitr_precision By = B[1];
    gitr_precision Bz = B[2];

    // Extract electric field components
    gitr_precision Ex = E[0];
    gitr_precision Ey = E[1];
    gitr_precision Ez = E[2];

    // Calculate q_over_m
    gitr_precision q_over_m = particlesPointer->charge[indx] * gitr_constants::e_charge / (particlesPointer->amu[indx] * gitr_constants::p_mass);

    // Half acceleration by electric field
    gitr_precision ux_minus = ux_i + 0.5 * q_over_m * Ex * particlesPointer->dt[indx];
    gitr_precision uy_minus = uy_i + 0.5 * q_over_m * Ey * particlesPointer->dt[indx];
    gitr_precision uz_minus = uz_i + 0.5 * q_over_m * Ez * particlesPointer->dt[indx];

    // Rotation by magnetic field
    gitr_precision t_x = q_over_m * Bx * 0.5 * particlesPointer->dt[indx];
    gitr_precision t_y = q_over_m * By * 0.5 * particlesPointer->dt[indx];
    gitr_precision t_z = q_over_m * Bz * 0.5 * particlesPointer->dt[indx];

    gitr_precision denom = 1 + t_x*t_x + t_y*t_y + t_z*t_z;
    gitr_precision s_x = 2*t_x / denom;
    gitr_precision s_y = 2*t_y / denom;
    gitr_precision s_z = 2*t_z / denom;

    gitr_precision ux_prime = ux_minus + uy_minus * t_z - uz_minus * t_y;
    gitr_precision uy_prime = uy_minus + uz_minus * t_x - ux_minus * t_z;
    gitr_precision uz_prime = uz_minus + ux_minus * t_y - uy_minus * t_x;

    gitr_precision ux_plus = ux_minus + uy_prime * s_z - uz_prime * s_y;
    gitr_precision uy_plus = uy_minus + uz_prime * s_x - ux_prime * s_z;
    gitr_precision uz_plus = uz_minus + ux_prime * s_y - uy_prime * s_x;

    // Half acceleration by electric field again
    gitr_precision ux_f = ux_plus + 0.5 * q_over_m * Ex * particlesPointer->dt[indx];
    gitr_precision uy_f = uy_plus + 0.5 * q_over_m * Ey * particlesPointer->dt[indx];
    gitr_precision uz_f = uz_plus + 0.5 * q_over_m * Ez * particlesPointer->dt[indx];

    // Updating particle positions
    gitr_precision x_f = x_i + ux_f * particlesPointer->dt[indx];
    gitr_precision y_f = y_i + uy_f * particlesPointer->dt[indx];
    gitr_precision z_f = z_i + uz_f * particlesPointer->dt[indx];

    // Finally, update particle properties directly
    particlesPointer->x[indx] = x_f;
    particlesPointer->y[indx] = y_f;
    particlesPointer->z[indx] = z_f;
    particlesPointer->vx[indx] = ux_f;
    particlesPointer->vy[indx] = uy_f;
    particlesPointer->vz[indx] = uz_f;

    return std::make_tuple(x_f, y_f, z_f, ux_f, uy_f, uz_f);
}
