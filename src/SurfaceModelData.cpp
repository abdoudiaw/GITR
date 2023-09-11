
// #include "Particles.h"
// #include "Surfaces.h"
// #include "array.h"
// #include "boris.h"
// #include "fieldLineTrace.h"
// #include "geometryCheck.h"
// #include "hashGeomSheath.h"
// #include "history.h"
// #include "interp2d.hpp"
// #include "interpolate.h"
// #include <cmath>
// #include "utils.h"
// #include <algorithm>
// #include <chrono>
// #include <cstdlib>
// #include <fstream>
// #include <iomanip>
// #include <iostream>
// #include <libconfig.h++>
// #include <netcdf.h>
// #include <random>
// #include <stdio.h>
// #include <stdlib.h>
// #include <vector>
// #include "flags.hpp"
// #include "surfaceModelData.h"


// struct surfaceReact
// {
//     /* data */
//     int nE_sputtRefCoeff, nA_sputtRefCoeff, nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOut, nE_sputtRefDistOutRef, nA_sputtRefDistOut, nDistE_surfaceModelRef, nDistE_surfaceModel, nDistA_surfaceModel;
//     std::string surfaceModelCfg, surfaceModelFile;
//     std::string input_path;
//     std::string cfgFile;
//     std::string surfaceModelCfg;
//     sim::Array<gitr_precision> E_sputtRefCoeff, A_sputtRefCoeff, Elog_sputtRefCoeff, energyDistGrid01, energyDistGrid01Ref, angleDistGrid01, spyl_surfaceModel, rfyl_surfaceModel, E_sputtRefDistIn, A_sputtRefDistIn, Elog_sputtRefDistIn, E_sputtRefDistOut, E_sputtRefDistOutRef, Aphi_sputtRefDistOut, Atheta_sputtRefDistOut, AphiDist_Y, AthetaDist_Y, EDist_Y, AphiDist_R, AthetaDist_R, EDist_R, AphiDist_CDF_Y, AthetaDist_CDF_Y, EDist_CDF_Y, AphiDist_CDF_R, AthetaDist_CDF_R, EDist_CDF_R, AphiDist_CDF_Y_regrid, AthetaDist_CDF_Y_regrid, EDist_CDF_Y_regrid, AphiDist_CDF_R_regrid, AthetaDist_CDF_R_regrid, EDist_CDF_R_regrid;cd ..
// };

// SurfaceModelData loadSurfaceModelData() {
//     SurfaceModelData data;
//     data.nE_sputtRefCoeff = 0;
//     data.nA_sputtRefCoeff = 0;
//     data.nE_sputtRefDistIn = 0;
//     data.nA_sputtRefDistIn = 0;
//     data.nE_sputtRefDistOut = 0;
//     data.nE_sputtRefDistOutRef = 0;
//     data.nA_sputtRefDistOut = 0;
//     data.nDistE_surfaceModelRef = 0;
//     data.nDistE_surfaceModel = 0;
//     data.nDistA_surfaceModel = 0;
//     data.surfaceModelCfg = "surfaceModel.cfg";
//     data.surfaceModelFile = "surfaceModel.nc";
//     data.input_path = "input/materialData/";
//     data.cfgFile = "input/surfaceModel.cfg";

//      const std::string input_path = "input/";
//     surfaceModelFile = "simpleSurfaceModel8ev.nc";
//     // config file is gitrInput.cfg
//     libconfig::Config cfg;
//     std::string cfgFile = input_path + "gitrInput.cfg";


//     getVariable(cfg, surfaceModelCfg + "fileString", surfaceModelFile);
//     nE_sputtRefCoeff = getDimFromFile(cfg, input_path + surfaceModelFile,
//                                         surfaceModelCfg, "nEsputtRefCoeffString");
//     nA_sputtRefCoeff = getDimFromFile(cfg, input_path + surfaceModelFile,
//                                         surfaceModelCfg, "nAsputtRefCoeffString");
//     nE_sputtRefDistIn =
//         getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                         "nEsputtRefDistInString");
//     nA_sputtRefDistIn =
//         getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                         "nAsputtRefDistInString");
//     nE_sputtRefDistOut =
//         getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                         "nEsputtRefDistOutString");
//     nE_sputtRefDistOutRef =
//         getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                         "nEsputtRefDistOutStringRef");
//     nA_sputtRefDistOut =
//         getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                         "nAsputtRefDistOutString");
//     nDistE_surfaceModel =
//         nE_sputtRefDistIn * nA_sputtRefDistIn * nE_sputtRefDistOut;
//     nDistE_surfaceModelRef =
//         nE_sputtRefDistIn * nA_sputtRefDistIn * nE_sputtRefDistOutRef;
//     nDistA_surfaceModel =
//         nE_sputtRefDistIn * nA_sputtRefDistIn * nA_sputtRefDistOut;
//     std::cout << " got dimensions of surface model " << std::endl;
    
//     sim::Array<gitr_precision> E_sputtRefCoeff(nE_sputtRefCoeff),
//         A_sputtRefCoeff(nA_sputtRefCoeff), Elog_sputtRefCoeff(nE_sputtRefCoeff),
//         energyDistGrid01(nE_sputtRefDistOut),
//         energyDistGrid01Ref(nE_sputtRefDistOutRef),
//         angleDistGrid01(nA_sputtRefDistOut),
//         spyl_surfaceModel(nE_sputtRefCoeff * nA_sputtRefCoeff),
//         rfyl_surfaceModel(nE_sputtRefCoeff * nA_sputtRefCoeff),
//         E_sputtRefDistIn(nE_sputtRefDistIn), A_sputtRefDistIn(nA_sputtRefDistIn),
//         Elog_sputtRefDistIn(nE_sputtRefDistIn),
//         E_sputtRefDistOut(nE_sputtRefDistOut),
//         E_sputtRefDistOutRef(nE_sputtRefDistOutRef),
//         Aphi_sputtRefDistOut(nA_sputtRefDistOut),
//         Atheta_sputtRefDistOut(nA_sputtRefDistOut),
//         AphiDist_Y(nDistA_surfaceModel), AthetaDist_Y(nDistA_surfaceModel),
//         EDist_Y(nDistE_surfaceModel), AphiDist_R(nDistA_surfaceModel),
//         AthetaDist_R(nDistA_surfaceModel), EDist_R(nDistE_surfaceModelRef),
//         AphiDist_CDF_Y(nDistA_surfaceModel),
//         AthetaDist_CDF_Y(nDistA_surfaceModel), EDist_CDF_Y(nDistE_surfaceModel),
//         AphiDist_CDF_R(nDistA_surfaceModel),
//         AthetaDist_CDF_R(nDistA_surfaceModel),
//         EDist_CDF_R(nDistE_surfaceModelRef),
//         AphiDist_CDF_Y_regrid(nDistA_surfaceModel),
//         AthetaDist_CDF_Y_regrid(nDistA_surfaceModel),
//         EDist_CDF_Y_regrid(nDistE_surfaceModel),
//         AphiDist_CDF_R_regrid(nDistA_surfaceModel),
//         AthetaDist_CDF_R_regrid(nDistA_surfaceModel),
//         EDist_CDF_R_regrid(nDistE_surfaceModelRef);


//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "E_sputtRefCoeff", E_sputtRefCoeff[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "A_sputtRefCoeff", A_sputtRefCoeff[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "E_sputtRefDistIn", E_sputtRefDistIn[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "A_sputtRefDistIn", A_sputtRefDistIn[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "E_sputtRefDistOut", E_sputtRefDistOut[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "E_sputtRefDistOutRef", E_sputtRefDistOutRef[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "Aphi_sputtRefDistOut", Aphi_sputtRefDistOut[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "Atheta_sputtRefDistOut", Atheta_sputtRefDistOut[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "sputtYldString", spyl_surfaceModel[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "reflYldString", rfyl_surfaceModel[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "EDist_Y", EDist_Y[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "AphiDist_Y", AphiDist_Y[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "AthetaDist_Y", AthetaDist_Y[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "EDist_R", EDist_R[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "AphiDist_R", AphiDist_R[0]);
//     getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
//                     "AthetaDist_R", AthetaDist_R[0]);

//     for (int i = 0; i < nE_sputtRefCoeff; i++) {
//         Elog_sputtRefCoeff[i] = log10(E_sputtRefCoeff[i]);
//         std::cout << " EsputtRefCoeff and Elog " << E_sputtRefCoeff[i] << " "
//                 << Elog_sputtRefCoeff[i] << std::endl;
//     }
//     for (int i = 0; i < nE_sputtRefDistIn; i++) {
//         Elog_sputtRefDistIn[i] = std::log10(E_sputtRefDistIn[i]);
//     }
//     for (int i = 0; i < nE_sputtRefDistOut; i++) {
//         energyDistGrid01[i] = i * 1.0 / nE_sputtRefDistOut;
//     }
//     for (int i = 0; i < nE_sputtRefDistOutRef; i++) {
//         energyDistGrid01Ref[i] = i * 1.0 / nE_sputtRefDistOutRef;
//     }
//     for (int i = 0; i < nA_sputtRefDistOut; i++) {
//         angleDistGrid01[i] = i * 1.0 / nA_sputtRefDistOut;
//     }
//     make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOut,
//                 EDist_Y.data(), EDist_CDF_Y.data());
//     make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 AphiDist_Y.data(), AphiDist_CDF_Y.data());
//     make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 AthetaDist_Y.data(), AthetaDist_CDF_Y.data());
//     make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOutRef,
//                 EDist_R.data(), EDist_CDF_R.data());
//     make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 AphiDist_R.data(), AphiDist_CDF_R.data());
//     make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 AthetaDist_R.data(), AthetaDist_CDF_R.data());
//     make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 AthetaDist_R.data(), AthetaDist_CDF_R.data());
//     regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 angleDistGrid01.data(), nA_sputtRefDistOut,
//                 Aphi_sputtRefDistOut[nA_sputtRefDistOut - 1],
//                 AphiDist_CDF_Y.data(), AphiDist_CDF_Y_regrid.data());
//     regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 angleDistGrid01.data(), nA_sputtRefDistOut,
//                 Atheta_sputtRefDistOut[nA_sputtRefDistOut - 1],
//                 AthetaDist_CDF_Y.data(), AthetaDist_CDF_Y_regrid.data());
//     regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOut,
//                 energyDistGrid01.data(), nE_sputtRefDistOut,
//                 E_sputtRefDistOut[nE_sputtRefDistOut - 1], EDist_CDF_Y.data(),
//                 EDist_CDF_Y_regrid.data());
//     regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 angleDistGrid01.data(), nA_sputtRefDistOut,
//                 Aphi_sputtRefDistOut[nA_sputtRefDistOut - 1],
//                 AphiDist_CDF_R.data(), AphiDist_CDF_R_regrid.data());
//     regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
//                 angleDistGrid01.data(), nA_sputtRefDistOut,
//                 Atheta_sputtRefDistOut[nA_sputtRefDistOut - 1],
//                 AthetaDist_CDF_R.data(), AthetaDist_CDF_R_regrid.data());
//     regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOutRef,
//                 energyDistGrid01Ref.data(), nE_sputtRefDistOutRef,
//                 E_sputtRefDistOutRef[nE_sputtRefDistOutRef - 1],
//                 EDist_CDF_R.data(), EDist_CDF_R_regrid.data());
    
//     gitr_precision spylAInterpVal = interp3d(
//         0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
//         nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
//         Elog_sputtRefDistIn.data(), AphiDist_CDF_Y_regrid.data());
//     gitr_precision spylAthetaInterpVal = interp3d(
//         0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
//         nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
//         Elog_sputtRefDistIn.data(), AthetaDist_CDF_Y_regrid.data());
//     gitr_precision sputEInterpVal = interp3d(
//         0.44, 63.0, std::log10(10.0), nE_sputtRefDistOut, nA_sputtRefDistIn,
//         nE_sputtRefDistIn, energyDistGrid01.data(), A_sputtRefDistIn.data(),
//         Elog_sputtRefDistIn.data(), EDist_CDF_Y_regrid.data());
//     gitr_precision rfylAInterpVal = interp3d(
//         0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
//         nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
//         Elog_sputtRefDistIn.data(), AphiDist_CDF_R_regrid.data());
//     gitr_precision rfylAthetaInterpVal = interp3d(
//         0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
//         nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
//         Elog_sputtRefDistIn.data(), AthetaDist_CDF_R_regrid.data());
//     gitr_precision rflEInterpVal = interp3d(
//         0.44, 63.0, std::log10(10.0), nE_sputtRefDistOut, nA_sputtRefDistIn,
//         nE_sputtRefDistIn, energyDistGrid01.data(), A_sputtRefDistIn.data(),
//         Elog_sputtRefDistIn.data(), EDist_CDF_R_regrid.data());

//     std::cout << "Finished surface model import sputtering" << spylAInterpVal
//                 << " " << spylAthetaInterpVal << " " << sputEInterpVal
//                 << std::endl;
//     std::cout << "Finished surface model import reflection" << rfylAInterpVal
//                 << " " << rfylAthetaInterpVal << " " << rflEInterpVal
//                 << std::endl;

//     // Update the surface model data
//     data.nE_sputtRefCoeff = nE_sputtRefCoeff;
//     data.nA_sputtRefCoeff = nA_sputtRefCoeff;
//     data.nE_sputtRefDistIn = nE_sputtRefDistIn;
//     data.nA_sputtRefDistIn = nA_sputtRefDistIn;
//     data.nE_sputtRefDistOut = nE_sputtRefDistOut;
//     data.nE_sputtRefDistOutRef = nE_sputtRefDistOutRef;
//     data.nA_sputtRefDistOut = nA_sputtRefDistOut;
//     data.nDistE_surfaceModelRef = nDistE_surfaceModelRef;
//     data.nDistE_surfaceModel = nDistE_surfaceModel;
//     data.nDistA_surfaceModel = nDistA_surfaceModel;
//     data.surfaceModelCfg = surfaceModelCfg;
//     data.surfaceModelFile = surfaceModelFile;
//     data.input_path = input_path;
//     data.cfgFile = cfgFile;
//     data.surfaceModelCfg = surfaceModelCfg;
//     data.E_sputtRefCoeff = E_sputtRefCoeff;
//     data.A_sputtRefCoeff = A_sputtRefCoeff;
//     data.Elog_sputtRefCoeff = Elog_sputtRefCoeff;
//     data.energyDistGrid01 = energyDistGrid01;
//     data.energyDistGrid01Ref = energyDistGrid01Ref;
//     data.angleDistGrid01 = angleDistGrid01;


// }