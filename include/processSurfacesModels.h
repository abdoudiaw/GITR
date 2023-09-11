// //------------------------------------------------------------------------------
// // GITR: processSurfacesModels.h
// //------------------------------------------------------------------------------
// //
// // Contributors:
// //     - GITR Community
// //
// // Last Modified:
// //     - August 2023 by Diaw
// //
// // Description:
// //     Processes surfaces reactions models data for reflection and sputtering and recombination. Requires 
// //     a file for each incident material and target material. The file should be named in the format: surfaceReactions_{incidentMaterial Z}_{targetMaterial Z}.nc
// //     and be located in the input/materialData     directory.
// // Note:
// //     This file is a component of the GITR codebase.
// //
// //------------------------------------------------------------------------------

// #include <netcdf>


// // Constants
// const gitr_precision SURFACE_BUFFER = 1.0e-6;

// // Cache for surface reaction data
// std::map<std::pair<std::string, std::string>, 
//          std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>> cache;

// inline std::tuple<size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> processSurfaceData(const std::string& incidentMaterial, const std::string& targetMaterial)
// {

//         auto key = std::make_pair(incidentMaterial, targetMaterial);
//         if (cache.find(key) != cache.end()) {
//             return cache[key];  // Return cached result
//         }

//     const std::string inputPath = "input/materialData/";
//     const std::string ratesFileName = "surfaceReactions_" + incidentMaterial + "_on_" + targetMaterial + ".nc";

//     // Open the netCDF file
//     netCDF::NcFile data(inputPath + ratesFileName, netCDF::NcFile::read);

//     printf("Loading surface reaction data from %s\n", ratesFileName.c_str());

//     // Get the variable objects
//     auto EVar = data.getVar("E");
//     auto AVar = data.getVar("A");
//     auto gridEnergyVar = data.getVar("nE");
//     auto gridAngleVar = data.getVar("nA");
//     auto sputterYieldVar = data.getVar("spyld");
//     auto reflectionYieldVar = data.getVar("rfyld");
//     auto cosXDistVar = data.getVar("cosXDist");
//     auto cosYDistVar = data.getVar("cosYDist");
//     auto cosXDistRefVar = data.getVar("cosXDistRef");
//     auto cosYDistRefVar = data.getVar("cosYDistRef");
//     auto energyDistVar = data.getVar("energyDist");
//     auto energyDistRefVar = data.getVar("energyDistRef");
//     auto eDistEgridVar = data.getVar("eDistEgrid");
//     auto eDistEgridRefVar = data.getVar("eDistEgridRef");
//     auto phiGridVar = data.getVar("phiGrid");
//     auto thetaGridVar = data.getVar("thetaGrid");

//     // Get the dimensions of the variables
//     auto energyDims = gridEnergyVar.getDims();
//     auto angleDims = gridAngleVar.getDims();

//     // Get the number of dimensions
//     size_t nEnergy = energyDims[0].getSize();
//     size_t nAngle = angleDims[0].getSize();

//     // Resize the arrays
//     std::vector<double> gridEnergy(nEnergy);
//     std::vector<double> gridAngle(nAngle);
//     std::vector<double> sputterYield(nEnergy * nAngle);
//     std::vector<double> reflectionYield(nEnergy * nAngle);

//     // Read the variable data into the arrays
//     gridEnergyVar.getVar(gridEnergy.data());
//     gridAngleVar.getVar(gridAngle.data());
//     sputterYieldVar.getVar(sputterYield.data());
//     reflectionYieldVar.getVar(reflectionYield.data());

//     // Close the netCDF file
//     data.close();

//     auto result = std::make_tuple(nEnergy, nAngle, gridEnergy, gridAngle, sputterYield, reflectionYield);
//     cache[key] = result;  // Store result in cache
//     return result;

// }
