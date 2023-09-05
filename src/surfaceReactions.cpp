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

// Constants
const gitr_precision SURFACE_BUFFER = 1.0e-6;

// Cache for surface reaction data
std::map<std::pair<std::string, std::string>, 
         std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>> cache;

std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> processSurfaceData(const std::string& incidentMaterial, const std::string& targetMaterial)
{
        auto key = std::make_pair(incidentMaterial, targetMaterial);
        if (cache.find(key) != cache.end()) {
            return cache[key];  // Return cached result
        }

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

    auto result = std::make_tuple(nEnergy, nAngle, gridEnergy, gridAngle, sputterYield, reflectionYield);
    cache[key] = result;  // Store result in cache
    return result;

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
     
        gitr_precision totalYR = 0.0;
        gitr_precision newWeight = 0.0;
        gitr_precision eInterpVal = 0.0;
        gitr_precision aInterpVal = 0.0;

        int wallHit = particles->surfaceHit[indx];
        int surfaceHit = boundaryVector[wallHit].surfaceNumber;
        int surface = boundaryVector[wallHit].surface;
        std::string material = materialData[boundaryVector[wallHit].Z].name;  
        std::string incidentName = particles->materialName[indx];
        gitr_precision surfaceBindingEnergy = materialData[boundaryVector[wallHit].Z].surfaceBindingEnergy;

        gitr_precision weight = particles->weight[indx];
        gitr_precision vx = particles->vx[indx];
        gitr_precision vy = particles->vy[indx];
        gitr_precision vz = particles->vz[indx];

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

            // Calculate the impact parameters
            calculateImpactParameters(thetaImpact, particleTrackVector, surfaceNormalVector);

            std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> surfaceData = processSurfaceData(particles->materialName[indx], material);

            size_t nEnergy = std::get<0>(surfaceData);
            size_t nAngle = std::get<1>(surfaceData);
            const std::vector<double>& gridEnergy = std::get<2>(surfaceData);
            const std::vector<double>& gridAngle = std::get<3>(surfaceData);
            const std::vector<double>& sputterYield = std::get<4>(surfaceData);
            const std::vector<double>& reflectionYield = std::get<5>(surfaceData);

            // Get the reflection probability
            gitr_precision reflectionProb = abs(interp2d(thetaImpact, std::log10(E0_for_surface_model), nAngle, nEnergy, gridAngle, gridEnergy, reflectionYield));

            // Get the sputtering probability
            gitr_precision sputteringProb = abs(interp2d(thetaImpact, std::log10(E0_for_surface_model), nAngle, nEnergy, gridAngle, gridEnergy, sputterYield));
            totalYR = reflectionProb + sputteringProb;

            // Get the maximum value in the energy grid
            gitr_precision Emax = gridEnergy[nEnergy - 1];
            auto sample = sampleThompsonDistribution(Emax, surfaceBindingEnergy);
            aInterpVal = sample.first;
            eInterpVal = sample.second;

            // Process the particle after impact
            processParticleAfterImpact(indx, particles, wallHit, newWeight, sputteringProb, totalYR, eInterpVal, aInterpVal, weight);
                    
        }
}

CUDA_CALLABLE_MEMBER_DEVICE
void calculateImpactParameters(gitr_precision& thetaImpact, const gitr_precision particleTrackVector[], const gitr_precision surfaceNormalVector[])
{
    gitr_precision partDotNormal = particleTrackVector[0] * surfaceNormalVector[0] + particleTrackVector[1] * surfaceNormalVector[1] + particleTrackVector[2] * surfaceNormalVector[2];
    thetaImpact = std::acos(partDotNormal);
    
    if (thetaImpact > M_PI / 2.0) {
        thetaImpact = std::abs(thetaImpact - M_PI);
    }
    thetaImpact = thetaImpact * 180.0 / M_PI;
    if (thetaImpact < 0.0) thetaImpact = 0.0;
}

CUDA_CALLABLE_MEMBER_DEVICE
void processParticleAfterImpact(int indx, Particles* particles,
    const int wallHit, gitr_precision& newWeight, const gitr_precision sputteringProb, const gitr_precision totalYR, gitr_precision& eInterpVal, gitr_precision& aInterpVal, gitr_precision weight)
{
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> distrib(0.0, 1.0);  // Define the range

    gitr_precision r7 = distrib(gen);

    if (totalYR > 0.0) {
        if (r7 > sputteringProb) {
            printf("Reflection\n");
            newWeight = weight * totalYR;
            // call reflect function
            // Reflect(indx, wallHit, newWeight, eInterpVal, aInterpVal, r7, boundaryVector, particles, flags);
            // Reflect (const int& indx, 
            //     const gitr_precision& newWeight,
            //     const gitr_precision& eInterpVal,
            //     const gitr_precision& aInterpVal,
            //     const gitr_precision& r10,
            //     Boundary* boundaryVector,
            //     Particles* particles,
            //     Flags* flags )


            //         particles->weight[indx] = newWeight;
        //     gitr_precision v[3];
        //     gitr_precision amu = particles->amu[indx];
        //     v[0] = particles->vx[indx];
        //     v[1] = particles->vy[indx];
        //     v[2] = particles->vz[indx];

        // // gitr_precision surfaceNormalVector[3] = {0.0};
        // boundaryVector[particles->surfaceHit[indx]].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

        // boundaryVector[particles->surfaceHit[indx]].transformToSurface(v, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

        // particles->vx[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vx;
        // particles->vy[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vy;
        // particles->vz[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vz;

        // particles->xprevious[indx] = particles->x[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[0] * SURFACE_BUFFER;
        // particles->yprevious[indx] = particles->y[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[1] * SURFACE_BUFFER;
        // particles->zprevious[indx] = particles->z[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[2] * SURFACE_BUFFER;
    

            // return \\\
        } else {
            printf("Sputtering\n");
            newWeight = weight * totalYR;
            if (sputteringProb == 0.0) newWeight = 0.0;
            // call sputter function
            // Sputter(wallHit, newWeight, eInterpVal, aInterpVal, r7, boundaryVector, particles, flags);
        }
    } else {
        newWeight = 0.0;
    }
}

// CUDA_CALLABLE_MEMBER_DEVICE
// void Sputter(const int& indx, 
//                 const gitr_precision& newWeight,
//                 const gitr_precision& eInterpVal,
//                 const gitr_precision& aInterpVal,
//                 const gitr_precision& r10,
//                 Boundary * boundaryVector,
//                 Particles* particles,
//                 Flags* flags ) {

//         gitr_precision amu = particles->amu[indx];
//         gitr_precision V0 = std::sqrt(2 * eInterpVal * gitr_constants::e_charge / (amu * gitr_constants::p_mass));
//         // set charge to zero
//         particles->charge[indx] = 0.0;

//         gitr_precision vSampled[3];
//         vSampled[0] = V0 * std::sin(aInterpVal * M_PI / 180) * std::cos(2.0 * M_PI * r10);
//         vSampled[1] = V0 * std::sin(aInterpVal * M_PI / 180) * std::sin(2.0 * M_PI * r10);
//         vSampled[2] = V0 * std::cos(aInterpVal * M_PI / 180);

//         gitr_precision surfaceNormalVector[3] = {0.0};
//         boundaryVector[particles->surfaceHit[indx]].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

//         boundaryVector[particles->surfaceHit[indx]].transformToSurface(vSampled, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

//         particles->vx[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vSampled[0];
//         particles->vy[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vSampled[1];
//         particles->vz[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vSampled[2];

//         particles->xprevious[indx] = particles->x[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[0] * SURFACE_BUFFER;
//         particles->yprevious[indx] = particles->y[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[1] * SURFACE_BUFFER;
//         particles->zprevious[indx] = particles->z[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[2] * SURFACE_BUFFER;
    
// }

// CUDA_CALLABLE_MEMBER_DEVICE
// void Reflect (const int& indx, 
//                 const gitr_precision& newWeight,
//                 const gitr_precision& eInterpVal,
//                 const gitr_precision& aInterpVal,
//                 const gitr_precision& r10,
//                 Boundary* boundaryVector,
//                 Particles* particles,
//                 Flags* flags ) {

//         particles->weight[indx] = newWeight;
//         gitr_precision amu = particles->amu[indx];
//         gitr_precision vx = particles->vx[indx];
//         gitr_precision vy = particles->vy[indx];
//         gitr_precision vz = particles->vz[indx];

//         gitr_precision v[3];
//         v[0] = vx;
//         v[1] = vy;
//         v[2] = vz;

//         gitr_precision surfaceNormalVector[3] = {0.0};
//         boundaryVector[particles->surfaceHit[indx]].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

//         boundaryVector[particles->surfaceHit[indx]].transformToSurface(v, particles->y[indx], particles->x[indx], flags->USE3DTETGEOM);

//         particles->vx[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vx;
//         particles->vy[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vy;
//         particles->vz[indx] = -static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * vz;

//         particles->xprevious[indx] = particles->x[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[0] * SURFACE_BUFFER;
//         particles->yprevious[indx] = particles->y[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[1] * SURFACE_BUFFER;
//         particles->zprevious[indx] = particles->z[indx] - static_cast<gitr_precision>(boundaryVector[particles->surfaceHit[indx]].inDir) * surfaceNormalVector[2] * SURFACE_BUFFER;
    
// }

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