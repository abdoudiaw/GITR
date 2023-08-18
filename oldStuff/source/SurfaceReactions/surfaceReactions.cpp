#include <iostream>
#include <tuple>
#include <vector>

#include <libconfig.h++>
#include <netcdf>
#include "thompson.h"
#include "surfaceReactions.h"
#include "constants.h"

enum ElementaryProcess
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
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
 int _nLines, Boundary* _boundaryVector, Surfaces* _surfaces, Flags* _flags) :
    particles(_particles), nLines(_nLines), boundaryVector(_boundaryVector), surfaces(_surfaces), state(_state), flags(_flags) { }

CUDA_CALLABLE_MEMBER_DEVICE
void surfaceReactions::operator()(std::size_t indx) const {
    if (particles->hitWall[indx] == 1.0) {
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
        gitr_precision eInterpVal = 0.0;
        gitr_precision aInterpVal = 0.0;
        gitr_precision weight = particles->weight[indx];
        gitr_precision vx = particles->vx[indx];
        gitr_precision vy = particles->vy[indx];
        gitr_precision vz = particles->vz[indx];
        gitr_precision SURFACE_BUFFER = 1.0e-4;

        if (flags->USE_SURFACEMODEL) {
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
            gitr_precision partDotNormal = vectorDotProduct(particleTrackVector, surfaceNormalVector);
            thetaImpact = std::acos(partDotNormal);
            if (thetaImpact > M_PI / 2.0) {
                thetaImpact = std::abs(thetaImpact - M_PI);
            }
            thetaImpact = thetaImpact * 180.0 / M_PI;
            if (thetaImpact < 0.0) thetaImpact = 0.0;
            // Get the data from the cache
            std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> surfaceData = processSurfaceData(particles->materialName[indx], std::string(1, surfaces->materialName[surfaceHit]));

            size_t nEnergy = std::get<0>(surfaceData);
            size_t nAngle = std::get<1>(surfaceData);
            const std::vector<double>& gridEnergy = std::get<2>(surfaceData);
            const std::vector<double>& gridAngle = std::get<3>(surfaceData);
            const std::vector<double>& sputterYield = std::get<4>(surfaceData);
            const std::vector<double>& reflectionYield = std::get<5>(surfaceData);



            // Get the reflection probability
            gitr_precision reflectionProb = interp2d(thetaImpact, std::log10(E0_for_surface_model), nAngle, nEnergy, gridAngle, gridEnergy, reflectionYield);
            // Get the sputtering probability
            gitr_precision sputteringProb = interp2d(thetaImpact, std::log10(E0_for_surface_model), nAngle, nEnergy, gridAngle, gridEnergy, sputterYield);
            totalYR = reflectionProb + sputteringProb;

            Y0 = sputteringProb;
            R0 = reflectionProb;


            // Get the maximum value in the energy grid
            gitr_precision Emax = gridEnergy[nEnergy - 1];
            auto sample = sampleThompsonDistribution(Emax, surfaces->bindingEnergy);
            aInterpVal = sample.first;
            eInterpVal = sample.second;
        
#ifdef __CUDACC__
            gitr_precision r7 = curand_uniform(&state[indx]);
            gitr_precision r8 = curand_uniform(&state[indx]);
            gitr_precision r9 = curand_uniform(&state[indx]);
            gitr_precision r10 = curand_uniform(&state[indx]);
#else
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            gitr_precision r7 = dist(state[indx]);
            gitr_precision r8 = dist(state[indx]);
            gitr_precision r9 = dist(state[indx]);
            gitr_precision r10 = dist(state[indx]);
#endif

            if (totalYR > 0.0) {
                if (dist(state[indx]) > sputteringProb) {
                    newWeight = weight * (totalYR);
                    if (surface > 0) {
#if USE_CUDA > 0
                        atomicAdd1(&surfaces->grossDeposition[surfaceHit], weight * (1.0 - R0));
                        atomicAdd1(&surfaces->grossErosion[surfaceHit], weight * Y0);
#else
                        surfaces->grossDeposition[surfaceHit] += (weight * (1.0 - R0));
                        surfaces->grossErosion[surfaceHit] += (weight * Y0);
#endif
                    }
                } else {
                    newWeight = weight * totalYR;
                    if (sputteringProb == 0.0) newWeight = 0.0;
                    if (surface > 0) {
#if USE_CUDA > 0
                        atomicAdd1(&surfaces->grossDeposition[surfaceHit], weight * (1.0 - R0));
                        atomicAdd1(&surfaces->grossErosion[surfaceHit], weight * Y0);
                        atomicAdd1(&surfaces->aveSputtYld[surfaceHit], Y0);
                        if (weight > 0.0) {
                            atomicAdd1(&surfaces->sputtYldCount[surfaceHit], 1.0);
                        }
#else
                        surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit] + weight * (1.0 - R0);
                        surfaces->grossErosion[surfaceHit] = surfaces->grossErosion[surfaceHit] + weight * Y0;
                        surfaces->aveSputtYld[surfaceHit] = surfaces->aveSputtYld[surfaceHit] + Y0;
                        surfaces->sputtYldCount[surfaceHit] = surfaces->sputtYldCount[surfaceHit] + 1;
#endif
                    }
                }
            } else {
                newWeight = 0.0;
                particles->hitWall[indx] = 2.0;
                if (surface > 0) {
#if USE_CUDA > 0
                    atomicAdd1(&surfaces->grossDeposition[surfaceHit], weight);
#else
                    surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit] + weight;
#endif
                }
            }

            if (eInterpVal <= 0.0) {
                newWeight = 0.0;
                particles->hitWall[indx] = 2.0;
                if (surface > 0) {
#if USE_CUDA > 0
                    atomicAdd1(&surfaces->grossDeposition[surfaceHit], weight * R0);
                    atomicAdd1(&surfaces->grossDeposition[surfaceHit], -weight * Y0);
#else
                    surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit] + weight;
#endif
                }
            }

            if (surface) {
#if USE_CUDA > 0
                atomicAdd1(&surfaces->sumWeightStrike[surfaceHit], weight);
                atomicAdd1(&surfaces->sumParticlesStrike[surfaceHit], 1.0);
#else
                surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit] + weight;
#endif
            }

            if (boundaryVector[wallHit].Z > 0.0 && newWeight > 0.0) {
                particles->weight[indx] = newWeight;
                particles->hitWall[indx] = 0.0;
                particles->charge[indx] = 0.0;
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
            } else {
                particles->hitWall[indx] = 2.0;
            }
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
