#include <iostream>
#include <fstream>
#include <string>
#include <netcdf.h>
#include <vector>


template<typename T>
std::vector<T> getVarData(const netCDF::NcVar& var) {
    std::vector<T> data(var.getDim(0).getSize());
    var.getVar(&data[0]);
    return data;
}

std::vector<std::string> getVarData(netCDF::NcFile& file, const std::string& varName) {
    int varId;
    nc_inq_varid(file.getId(), varName.c_str(), &varId);
    nc_type varType;
    nc_inq_vartype(file.getId(), varId, &varType);
    if (varType != NC_STRING) {
        throw std::runtime_error("The variable is not of string type.");
    }
    size_t varLen;
    nc_inq_dimlen(file.getId(), 0, &varLen); // Assuming 0 is the dimension ID for 'nP'
    std::vector<char*> buffer(varLen);
    nc_get_var_string(file.getId(), varId, buffer.data());
    std::vector<std::string> result(varLen);
    for (size_t i = 0; i < varLen; ++i) {
        result[i] = buffer[i];
        free(buffer[i]); // Free memory allocated by nc_get_var_string
    }
    return result;
}


struct ParticleData {
    std::vector<gitr_precision> impurity_id, xpfile, ypfile, zpfile,
        vxpfile, vypfile, vzpfile, charge, mass, ionizationState;
    std::vector<std::string> materialName;
};

ParticleData readParticleData(const std::string &ncParticleSourceFile) {
    // Read in particle source file
    ParticleData data;
    int nPfile = 0;

    std::cout << "Starting psourcefile import " << std::endl;
    
    // Return this in event of a problem.
    static const int NC_ERR = 2;

    try {
        netCDF::NcFile ncp0("input/" + ncParticleSourceFile, netCDF::NcFile::read);
    }
    catch (netCDF::exceptions::NcException &e)
    {
        e.what();
        std::cout << "FAILURE*************************************" << std::endl;
        throw std::runtime_error("Failed to open NetCDF file.");
    }
        
    std::cout << "finished NcFile ncp0 starting ncp" << std::endl;
    netCDF::NcFile ncp("input/" + ncParticleSourceFile, netCDF::NcFile::read);
    std::cout << "getting dim nP" << std::endl;

    // char material[100];
    netCDF::NcDim ps_nP(ncp.getDim("nP"));
    nPfile = ps_nP.getSize();

    // Initialize data vectors with correct size
    data.impurity_id.resize(nPfile); data.xpfile.resize(nPfile); data.ypfile.resize(nPfile);
    data.zpfile.resize(nPfile); data.vxpfile.resize(nPfile); data.vypfile.resize(nPfile);
    data.vzpfile.resize(nPfile); data.charge.resize(nPfile); data.mass.resize(nPfile);
    data.ionizationState.resize(nPfile); data.materialName.resize(nPfile);

    // Get variable references
    netCDF::NcVar ncp_impurity_id(ncp.getVar("impurity_id"));
    netCDF::NcVar ncp_x(ncp.getVar("x"));
    netCDF::NcVar ncp_y(ncp.getVar("y"));
    netCDF::NcVar ncp_z(ncp.getVar("z"));
    netCDF::NcVar ncp_vx(ncp.getVar("vx"));
    netCDF::NcVar ncp_vy(ncp.getVar("vy"));
    netCDF::NcVar ncp_vz(ncp.getVar("vz"));
    netCDF::NcVar ncp_charge(ncp.getVar("charge"));
    netCDF::NcVar ncp_mass(ncp.getVar("mass"));
    netCDF::NcVar ncp_IonizationState(ncp.getVar("IonizationState"));
    netCDF::NcVar ncp_materialName(ncp.getVar("materialName"));

    // Get variable data
    data.impurity_id = getVarData<gitr_precision>(ncp_impurity_id);
    data.xpfile = getVarData<gitr_precision>(ncp_x);
    data.ypfile = getVarData<gitr_precision>(ncp_y);
    data.zpfile = getVarData<gitr_precision>(ncp_z);
    data.vxpfile = getVarData<gitr_precision>(ncp_vx);
    data.vypfile = getVarData<gitr_precision>(ncp_vy);
    data.vzpfile = getVarData<gitr_precision>(ncp_vz);
    data.charge = getVarData<gitr_precision>(ncp_charge);
    data.mass = getVarData<gitr_precision>(ncp_mass);
    data.ionizationState = getVarData<gitr_precision>(ncp_IonizationState);
    data.materialName = getVarData(ncp, "materialName");
    ncp.close();
    return data;

    }

    void initializeParticleArray(const ParticleData& particleData, Particles* particleArray,
    sim::Array<gitr_precision>& px, sim::Array<gitr_precision>& py, sim::Array<gitr_precision>& pz,
    sim::Array<gitr_precision>& pvx, sim::Array<gitr_precision>& pvy, sim::Array<gitr_precision>& pvz,
    sim::Array<gitr_precision>& pZ, sim::Array<gitr_precision>& pamu, sim::Array<gitr_precision>& pcharge,
    gitr_precision dt,  sim::Array<std::string >& pmaterialName )
    {
    int nP = particleData.materialName.size();
    long nParticles = nP;

    auto& xpfile = particleData.xpfile;
    auto& ypfile = particleData.ypfile;
    auto& zpfile = particleData.zpfile;
    auto& vxpfile = particleData.vxpfile;
    auto& vypfile = particleData.vypfile;
    auto& vzpfile = particleData.vzpfile;
    auto& charge = particleData.charge;
    auto& mass = particleData.mass;
    auto& ionizationState = particleData.ionizationState;
    auto& materialName = particleData.materialName;

    for (size_t i = 0; i < nParticles; ++i) {
        // Set particle data
        gitr_precision x = xpfile[i];
        gitr_precision y = ypfile[i];
        gitr_precision z = zpfile[i];
        gitr_precision vx = vxpfile[i];
        gitr_precision vy = vypfile[i];
        gitr_precision vz = vzpfile[i];
        gitr_precision particleCharge = charge[i];
        gitr_precision particleMass = mass[i];
        gitr_precision ionization = ionizationState[i];
        std::string material = materialName[i];

        particleArray->setParticleV(i, x, y, z, vx, vy, vz, ionization, particleMass, particleCharge, dt, material );
 
        px[i] = x;
        py[i] = y;
        pz[i] = z;
        pvx[i] = vx;
        pvy[i] = vy;
        pvz[i] = vz;
        pZ[i] = ionization;
        pamu[i] = particleMass;
        pcharge[i] = particleCharge;
        pmaterialName[i] = material;

    }

}


// Function to process plasma profiles
std::tuple<int, int, int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> processPlasmaBackground(const std::string &filename)
{
    std::string input_path = "input/";
//    std::string plasmaFile = "plasmaProfile.nc";
    std::string plasmaFile = filename;

    // Open the netCDF file
//    netCDF::NcFile data(filename, netCDF::NcFile::read);
    netCDF::NcFile data(input_path + plasmaFile, netCDF::NcFile::read);

    // Get the variable objects
    netCDF::NcVar gridTe_var = data.getVar("te");
    netCDF::NcVar gridNe_var = data.getVar("ne");
    netCDF::NcVar gridTi_var = data.getVar("ti");
    netCDF::NcVar gridNi_var = data.getVar("ni");
    netCDF::NcVar gridVx_var = data.getVar("vx");
    netCDF::NcVar gridVy_var = data.getVar("vy");
    netCDF::NcVar gridVz_var = data.getVar("vz");

    // Get the dimensions of the variables
    std::vector<netCDF::NcDim> temperatureDims = gridTe_var.getDims();
    std::vector<netCDF::NcDim> densityDims = gridNe_var.getDims();
    std::vector<netCDF::NcDim> velocityDimsX = gridVx_var.getDims();
    std::vector<netCDF::NcDim> velocityDimsY = gridVy_var.getDims();
    std::vector<netCDF::NcDim> velocityDimsZ = gridVz_var.getDims();

    // Get the number of temperatures, densities, and velocities
    size_t nTemperatures = temperatureDims[0].getSize();
    size_t nDensities = densityDims[0].getSize();
    size_t nVelocitiesX = velocityDimsX[0].getSize();
    size_t nVelocitiesY = velocityDimsY[0].getSize();
    size_t nVelocitiesZ = velocityDimsZ[0].getSize();

    // Resize the arrays
    std::vector<double> gridTemperatureElectrons(nTemperatures);
    std::vector<double> gridTemperatureIons(nTemperatures);
    std::vector<double> grideDensity(nDensities);
    std::vector<double> gridiDensity(nDensities);
    std::vector<double> gridVelocityX(nVelocitiesX);
    std::vector<double> gridVelocityY(nVelocitiesY);
    std::vector<double> gridVelocityZ(nVelocitiesZ);

    // Read the variable data into the arrays
    gridTe_var.getVar(gridTemperatureElectrons.data());
    gridTi_var.getVar(gridTemperatureIons.data());
    gridNe_var.getVar(grideDensity.data());
    gridNi_var.getVar(gridiDensity.data());
    gridVx_var.getVar(gridVelocityX.data());
    gridVy_var.getVar(gridVelocityY.data());
    gridVz_var.getVar(gridVelocityZ.data());

    // Close the netCDF file
    data.close();

    return std::make_tuple(nTemperatures, nDensities, nVelocitiesX, std::move(gridTemperatureElectrons), std::move(gridTemperatureIons), std::move(grideDensity), std::move(gridiDensity),
    std::move(gridVelocityX), std::move(gridVelocityY), std::move(gridVelocityZ));
}



