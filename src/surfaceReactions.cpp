#include <iostream>
#include <tuple>
#include <vector>

#include <libconfig.h++>
#include <netcdf>
#include "constants.h"
#include "thompson.h"
#include "surfaceReactions.h"
#include "targetNames.h"
#include <random>


enum surfacesReactions
{
    REFLECTION,
    SPUTTERING
};

std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> processSurfaceData(const std::string& incidentMaterial, const std::string& targetMaterial)
{
    const std::string inputPath = "input/materialData/";
    const std::string ratesFileName = "surfaceReactions_" + incidentMaterial + "_on_" + targetMaterial + ".nc";


    // Open the netCDF file
    netCDF::NcFile data(inputPath + ratesFileName, netCDF::NcFile::read);

    printf("Loading surface reaction data from %s\n", ratesFileName.c_str());

    // Get the variable objects
    auto gridEnergyVar = data.getVar("nE");
    auto gridAngleVar = data.getVar("nA");
    auto sputterYieldVar = data.getVar("spyld");
    auto reflectionYieldVar = data.getVar("rfyld");

    // Get the dimensions of the variables
    auto energyDims = gridEnergyVar.getDims();
    auto angleDims = gridAngleVar.getDims();

    // Get the number of dimensions
    size_t nEnergy = energyDims[0].getSize();
    size_t nAngle = angleDims[0].getSize();

    // Resize the arrays
    std::vector<double> gridEnergy(nEnergy);
    std::vector<double> gridAngle(nAngle);
    std::vector<double> sputterYield(nEnergy * nAngle);
    std::vector<double> reflectionYield(nEnergy * nAngle);

    // Read the variable data into the arrays
    gridEnergyVar.getVar(gridEnergy.data());
    gridAngleVar.getVar(gridAngle.data());
    sputterYieldVar.getVar(sputterYield.data());
    reflectionYieldVar.getVar(reflectionYield.data());

    // Close the netCDF file
    data.close();

    return std::make_tuple(nEnergy, nAngle, gridEnergy, gridAngle, sputterYield, reflectionYield);
}

surfaceReactions::surfaceReactions(Particles* _particles,
 int _nLines, int _nSurfaces,
 Boundary* _boundaryVector, Flags* _flags) :
    particles(_particles), nLines(_nLines), nSurfaces(_nSurfaces), boundaryVector(_boundaryVector),  flags(_flags) { }


CUDA_CALLABLE_MEMBER_DEVICE
void surfaceReactions::operator()(std::size_t indx) const {

    // if (particles->hitWall[indx] == 1.0) {
        gitr_precision thetaImpact = 0.0;
        gitr_precision particleTrackVector[3] = {0.0};
        gitr_precision surfaceNormalVector[3] = {0.0};
        gitr_precision vSampled[3] = {0.0};
        gitr_precision Y0 = 0.0;
        gitr_precision R0 = 0.0;
        gitr_precision totalYR = 0.0;
        gitr_precision newWeight = 0.0;
        int wallHit = particles->surfaceHit[indx];
        int surfaceHit = boundaryVector[wallHit].surfaceNumber;
        int surface = boundaryVector[wallHit].surface;
        std::string material = materialData[boundaryVector[wallHit].Z].name;  
        std::string incidentName = particles->materialName[indx];
        gitr_precision surfaceBindingEnergy = materialData[boundaryVector[wallHit].Z].surfaceBindingEnergy;

        gitr_precision eInterpVal = 0.0;
        gitr_precision aInterpVal = 0.0;
        
        gitr_precision weight = particles->weight[indx];

        gitr_precision vx = particles->vx[indx];
        gitr_precision vy = particles->vy[indx];
        gitr_precision vz = particles->vz[indx];

        gitr_precision SURFACE_BUFFER = 1.0e-6;

        // Need to call this only once for each material and cached the data
        if (boundaryVector[wallHit].Z > 0)   std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> surfaceData = processSurfaceData(incidentName, material);


        if (boundaryVector[wallHit].Z > 0) {
            particles->firstCollision[indx] = 1;
            particleTrackVector[0] = vx;
            particleTrackVector[1] = vy;
            particleTrackVector[2] = vz;
            gitr_precision velocityNorm = std::sqrt(vx * vx + vy * vy + vz * vz);
            gitr_precision E0_for_surface_model = 0.5 * particles->amu[indx] * gitr_constants::p_mass * (velocityNorm * velocityNorm) / gitr_constants::e_charge;

            boundaryVector[wallHit].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

            particleTrackVector[0] = vx / velocityNorm;
            particleTrackVector[1] = vy / velocityNorm;
            particleTrackVector[2] = vz / velocityNorm;

            gitr_precision partDotNormal = particleTrackVector[0] * surfaceNormalVector[0] + particleTrackVector[1] * surfaceNormalVector[1] + particleTrackVector[2] * surfaceNormalVector[2];
            thetaImpact = std::acos(partDotNormal);

            if (thetaImpact > M_PI / 2.0) {
                thetaImpact = std::abs(thetaImpact - M_PI);
            }
            thetaImpact = thetaImpact * 180.0 / M_PI;
            if (thetaImpact < 0.0) thetaImpact = 0.0;
            // Get the data from the cache
    
            size_t nEnergy = 100; //std::get<0>(surfaceData);
            size_t nAngle = 10; //std::get<1>(surfaceData);
            // const std::vector<double>& gridEnergy = std::get<2>(surfaceData);
            // const std::vector<double>& gridAngle = std::get<3>(surfaceData);
            // const std::vector<double>& sputterYield = std::get<4>(surfaceData);
            // const std::vector<double>& reflectionYield = std::get<5>(surfaceData);

            // Get the reflection probability
            gitr_precision reflectionProb = 0.5; //interp2d(thetaImpact, std::log10(E0_for_surface_model), nAngle, nEnergy, gridAngle, gridEnergy, reflectionYield);
            // Get the sputtering probability
            gitr_precision sputteringProb = 0.5; //interp2d(thetaImpact, std::log10(E0_for_surface_model), nAngle, nEnergy, gridAngle, gridEnergy, sputterYield);
            totalYR = reflectionProb + sputteringProb;

            /// debugging case for testing
            Y0 = 0.5; 
            R0 = 0.5; 

            reflectionProb = R0 ;
            sputteringProb = Y0;
            totalYR = R0 + Y0; 

            // Get the maximum value in the energy grid
            // gitr_precision Emax = gridEnergy[nEnergy - 1];
            // auto sample = sampleThompsonDistribution(Emax, surfaceBindingEnergy);
            // aInterpVal = sample.first;
            // eInterpVal = sample.second;

            // set surface constructor
            auto surfaces = new Surfaces(nSurfaces, nEnergy, nAngle);
            // build grid for surface
            // surfaces->gridE = gridEnergy;
            // surfaces->gridA = gridAngle;

            std::random_device rd;  
            std::mt19937 gen(rd()); 
            std::uniform_real_distribution<> distrib(0.0, 1.0);  // Define the range

            gitr_precision r7 = distrib(gen);
            gitr_precision r8 = distrib(gen);
            gitr_precision r9 = distrib(gen);
            gitr_precision r10 = distrib(gen);

            if (totalYR > 0.0) {
                if (r7 > sputteringProb) {
                    newWeight = weight * (totalYR);
                } else {
                    newWeight = weight * totalYR;
                    if (sputteringProb == 0.0) newWeight = 0.0;
                }
            } else {
                newWeight = 0.0;
                particles->hitWall[indx] = 2.0;
            }
            // }
            if (boundaryVector[wallHit].Z > 0.0 && newWeight > 0.0) {
                particles->weight[indx] = newWeight;
                particles->hitWall[indx] = 0.0;
                particles->charge[indx] = 0.0;
                particles->Z[indx] = boundaryVector[wallHit].Z;
                particles->amu[indx] = boundaryVector[wallHit].amu;
                particles->materialName[indx] = material;
                gitr_precision V0 = std::sqrt(2 * eInterpVal * gitr_constants::e_charge / (particles->amu[indx] * gitr_constants::p_mass));

                vSampled[0] = V0 * std::sin(aInterpVal * M_PI / 180) * std::cos(2.0 * M_PI * r10);
                vSampled[1] = V0 * std::sin(aInterpVal * M_PI / 180) * std::sin(2.0 * M_PI * r10);
                vSampled[2] = V0 * std::cos(aInterpVal * M_PI / 180);

                boundaryVector[wallHit].transformToSurface(vSampled, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

                particles->vx[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[0];
                particles->vy[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[1];
                particles->vz[indx] = -static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * vSampled[2];
                
                particles->xprevious[indx] = particles->x[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[0] * SURFACE_BUFFER;
                particles->yprevious[indx] = particles->y[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[1] * SURFACE_BUFFER;
                particles->zprevious[indx] = particles->z[indx] - static_cast<gitr_precision>(boundaryVector[wallHit].inDir) * surfaceNormalVector[2] * SURFACE_BUFFER;
            } 
        }
}

 CUDA_CALLABLE_MEMBER
 gitr_precision interp2d(gitr_precision x, gitr_precision z, int nx, int nz,
                         const std::vector<double>& gridx, const std::vector<double>& gridz, const std::vector<double>& data)
 {
     gitr_precision fxz = 0.0;
     gitr_precision fx_z1 = 0.0;
     gitr_precision fx_z2 = 0.0;

     if (nx * nz == 1) {
         fxz = data[0];
     } else {
         gitr_precision dim1 = x;
         gitr_precision d_dim1 = gridx[1] - gridx[0];
         gitr_precision dz = gridz[1] - gridz[0];
         int i = std::floor((dim1 - gridx[0]) / d_dim1);
         int j = std::floor((z - gridz[0]) / dz);
         if (i < 0) i = 0;
         if (j < 0) j = 0;
         if (i >= nx - 1 && j >= nz - 1) {
             fxz = data[nx - 1 + (nz - 1) * nx];
         } else if (i >= nx - 1) {
             fx_z1 = data[nx - 1 + j * nx];
             fx_z2 = data[nx - 1 + (j + 1) * nx];
             fxz = ((gridz[j + 1] - z) * fx_z1 + (z - gridz[j]) * fx_z2) / dz;
         } else if (j >= nz - 1) {
             fx_z1 = data[i + (nz - 1) * nx];
             fx_z2 = data[i + (nz - 1) * nx];
             fxz = ((gridx[i + 1] - dim1) * fx_z1 + (dim1 - gridx[i]) * fx_z2) / d_dim1;
         } else {
             fx_z1 = ((gridx[i + 1] - dim1) * data[i + j * nx] + (dim1 - gridx[i]) * data[i + 1 + j * nx]) / d_dim1;
             fx_z2 = ((gridx[i + 1] - dim1) * data[i + (j + 1) * nx] + (dim1 - gridx[i]) * data[i + 1 + (j + 1) * nx]) / d_dim1;
             fxz = ((gridz[j + 1] - z) * fx_z1 + (z - gridz[j]) * fx_z2) / dz;
         }
     }

     return fxz;
 }

