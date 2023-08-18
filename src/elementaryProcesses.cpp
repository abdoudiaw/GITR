/* Copyright 2023
 *
 * This file is part of GITR.
 *
 * License: BSD-3-Clause-LBNL
 * Website: https://github.com/ORNL-FUSION/GITR/tree/dev
 * Modified by on 05-20-2023 15:17 by Abdou Diaw
 */

#include "elementaryProcesses.h"
#include "interpolater.h"
#include <libconfig.h++>
#include <netcdf>

enum ElementaryProcess
{
    IONIZATION,
    RECOMBINATION
};

#if __CUDA_ARCH__
CUDA_CALLABLE_MEMBER_DEVICE
double get_rand_double(curandState* state, int indx)
{
    curandState localState = state[indx];
    double r1 = curand_uniform_double(&localState);
    state[indx] = localState;
    return r1;
}
#endif

double get_rand_double(std::mt19937* state, int indx)
{
    std::uniform_real_distribution<double> dist(-1.0, 0.0);
    std::mt19937 this_state = state[indx];
    double r1 = -dist(this_state);
    state[indx] = this_state;
    return r1;
}


/* Function to process ADAS data for the elementary processes listed above */
std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>> process_rates(ElementaryProcess process, int charge)
{

    std::string filePrefix;
    if (process == IONIZATION)
    {
        filePrefix = "Ionization";
    }
    else if (process == RECOMBINATION)
    {
        filePrefix = "Recombination";
    }
    std::string input_path = "input/adasData/";
    std::string ratesFiles = "ADAS_Rates_" + std::to_string(charge) + ".nc";

    // Open the netCDF file
    netCDF::NcFile data(input_path + ratesFiles, netCDF::NcFile::read);

    // Get the variable objects
    netCDF::NcVar gridTemperature_var = data.getVar("gridTemperature_" + filePrefix);
    netCDF::NcVar gridDensity_var = data.getVar("gridDensity_" + filePrefix);
    netCDF::NcVar rateCoeff_var = data.getVar(filePrefix + "RateCoeff");

    // Get the dimensions of the variables
    std::vector<netCDF::NcDim> temperatureDims = gridTemperature_var.getDims();
    std::vector<netCDF::NcDim> densityDims = gridDensity_var.getDims();
    std::vector<netCDF::NcDim> chargeStateDims = rateCoeff_var.getDims();

    // Get the number of temperatures, densities, and charge states
    size_t nTemperatures = temperatureDims[0].getSize();
    size_t nDensities = densityDims[0].getSize();
    size_t nChargeStates = chargeStateDims[0].getSize();

    // Resize the arrays
    std::vector<double> gridTemperature(nTemperatures);
    std::vector<double> gridDensity(nDensities);
    std::vector<double> rateCoeff(nChargeStates * nTemperatures * nDensities);

    // Read the variable data into the arrays
    gridTemperature_var.getVar(gridTemperature.data());
    gridDensity_var.getVar(gridDensity.data());
    rateCoeff_var.getVar(rateCoeff.data());

    // Close the netCDF file
    data.close();

    return std::make_tuple(nTemperatures, nDensities, nChargeStates, gridTemperature, gridDensity, rateCoeff);
}

template<typename T>
elementaryProcesses<T>::elementaryProcesses(Particles* _particlesPointer, T* _state, int _nR_Dens, int _nZ_Dens, gitr_precision* _DensGridr, gitr_precision* _DensGridz, gitr_precision* _ne, int _nR_Temp, int _nZ_Temp, gitr_precision* _TempGridr, gitr_precision* _TempGridz, gitr_precision* _te, Flags* _flags)
    : particlesPointer(_particlesPointer), nR_Dens(_nR_Dens), nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz), ne(_ne), nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), TempGridz(_TempGridz), te(_te), state(_state),  flags (_flags)
{ }

template<typename T>
CUDA_CALLABLE_MEMBER_DEVICE
void elementaryProcesses<T>::operator()(std::size_t indx)
{
    if (flags->USE_IONIZATION)
    {
        gitr_precision dt = particlesPointer->dt[indx];

        // Get the data from the cache
        std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>> ionizationResult = process_rates(IONIZATION, particlesPointer->Z[indx]);
        size_t nTemperaturesIonize = std::get<0>(ionizationResult);
        size_t nDensitiesIonize = std::get<1>(ionizationResult);
        size_t nChargeStatesIonize = std::get<2>(ionizationResult);
        const std::vector<double>& gridTemperature_Ionization = std::get<3>(ionizationResult);
        const std::vector<double>& gridDensity_Ionization = std::get<4>(ionizationResult);
        const std::vector<double>& rateCoeff_Ionization = std::get<5>(ionizationResult);

        // Get local temperature and density
        gitr_precision tlocal = interp2dCombined(particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],
                                                 nR_Temp, nZ_Temp, TempGridr, TempGridz, te);
        gitr_precision nlocal = interp2dCombined(particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],
            nR_Dens, nZ_Dens, DensGridr, DensGridz, ne);

        // Get the ionization rate coefficient
        gitr_precision RClocal = rateCoeffInterp(particlesPointer->charge[indx], tlocal, nlocal, nTemperaturesIonize, nDensitiesIonize, gridTemperature_Ionization, gridDensity_Ionization, rateCoeff_Ionization);
        tion = 1.0 / (RClocal * nlocal);
        if (tlocal == 0.0 || nlocal == 0.0) tion = 1.0e12;

        gitr_precision P = exp(-dt / tion);
        gitr_precision P1 = 1.0 - P;
        gitr_precision r1 = get_rand_double(state, indx);

        if (particlesPointer->hitWall[indx] == 0.0)
        {
            if (r1 <= P1)
            {
                particlesPointer->charge[indx] = particlesPointer->charge[indx] + 1;
                std::cout << "Ionization" << std::endl;
            }
        }
        
    }

    if (flags->USE_RECOMBINATION)
    {
        gitr_precision dt = particlesPointer->dt[indx];

        // Get the data from the cache
        std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>> recombinationResult = process_rates(RECOMBINATION, particlesPointer->Z[indx]);
        size_t nTemperaturesRecomb = std::get<0>(recombinationResult);
        size_t nDensitiesRecomb = std::get<1>(recombinationResult);
        size_t nChargeStatesRecomb = std::get<2>(recombinationResult);
        const std::vector<double>& gridTemperature_Recombination = std::get<3>(recombinationResult);
        const std::vector<double>& gridDensity_Recombination = std::get<4>(recombinationResult);
        const std::vector<double>& rateCoeff_Recombination = std::get<5>(recombinationResult);

        // Get local temperature and density
        gitr_precision tlocal = interp2dCombined(particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],
            nR_Temp, nZ_Temp, TempGridr, TempGridz, te);
        gitr_precision nlocal = interp2dCombined(particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],
            nR_Dens, nZ_Dens, DensGridr, DensGridz, ne);

        // Get the recombination rate coefficient
        gitr_precision RClocal = rateCoeffInterp(particlesPointer->charge[indx], tlocal, nlocal, nTemperaturesRecomb, nDensitiesRecomb, gridTemperature_Recombination, gridDensity_Recombination, rateCoeff_Recombination);
        gitr_precision trecomb = 1.0 / (RClocal * nlocal);
        if (tlocal == 0.0 || nlocal == 0.0) trecomb = 1.0e12;


        gitr_precision P = exp(-dt/trecomb);
        gitr_precision P1 = 1.0-P;
        gitr_precision r1 = get_rand_double(state, indx);

        if(particlesPointer->hitWall[indx] == 0.0)
        {
          if(r1 <= P1)
          {
            particlesPointer->charge[indx] = particlesPointer->charge[indx]-1;
            std::cout << "Recombination" << std::endl;
          }
        }
          
    }
}

/* This function interpolates the rate coefficient from the rate coefficient grid */
gitr_precision rateCoeffInterp(int charge, gitr_precision te, gitr_precision ne, int nT, int nD, const std::vector<double>& rateGrid_Tempp, const std::vector<double>& rateGrid_Densp, const std::vector<double>& Ratesp)
{
    int indT = 0;
    int indN = 0;
    gitr_precision logT = std::log10(te);
    gitr_precision logn = std::log10(ne);
    gitr_precision d_T = rateGrid_Tempp[1] - rateGrid_Tempp[0];
    gitr_precision d_n = rateGrid_Densp[1] - rateGrid_Densp[0];
    indT = std::floor((logT - rateGrid_Tempp[0]) / d_T);
    indN = std::floor((logn - rateGrid_Densp[0]) / d_n);
    if (indT < 0 || indT > nT - 2) {
        indT = 0;
    }
    if (indN < 0 || indN > nD - 2) {
        indN = 0;
    }
    if (charge > 74 - 1) {
        charge = 0;
    }
    gitr_precision aT = std::pow(10.0, rateGrid_Tempp[indT + 1]) - te;
    gitr_precision bT = te - std::pow(10.0, rateGrid_Tempp[indT]);
    gitr_precision abT = aT + bT;
    gitr_precision aN = std::pow(10.0, rateGrid_Densp[indN + 1]) - ne;
    gitr_precision bN = ne - std::pow(10.0, rateGrid_Densp[indN]);
    gitr_precision abN = aN + bN;
    gitr_precision fx_z1 = (aN * std::pow(10.0, Ratesp[charge * nT * nD + indT * nD + indN]) +
        bN * std::pow(10.0, Ratesp[charge * nT * nD + indT * nD + indN + 1])) / abN;
    gitr_precision fx_z2 = (aN * std::pow(10.0, Ratesp[charge * nT * nD + (indT + 1) * nD + indN]) +
        bN * std::pow(10.0, Ratesp[charge * nT * nD + (indT + 1) * nD + indN + 1])) / abN;
    gitr_precision fxz = (aT * fx_z1 + bT * fx_z2) / abT;
    return fxz;
}

/* explicit template instantiations */
#if USE_CUDA
template struct elementaryProcesses<curandState>;
#else
template struct elementaryProcesses<std::mt19937>;
#endif

