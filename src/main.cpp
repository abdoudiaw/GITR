#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libconfig.h++>
#include <netcdf.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <set>

#include </opt/homebrew/opt/libomp/include/omp.h>
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include "array.h"
#include "boundaryInit.h"
#include "boundary.h"
#include "configInterface.h"
#include "curandInitialize.h"
#include "flags.h"
#include "backgroundCollisionsMCC.h"
#include "geometryCheck.h"
#include "interpolater.h"
#include "elementaryProcesses.h"
#include "particles.h"
#include "pusher.h"
#include "surfaceReactions.h"
#include "utils.h"
#include "dump.h"
#include "getParticleData.h"

#ifdef __CUDACC__
#include <curand.h>
#include <curand_kernel.h>
#endif
#include "CLI/CLI.hpp"
using namespace netCDF;
#if USE_DOUBLE
typedef double gitr_precision;
netCDF::NcType netcdf_precision = netCDF::ncDouble;
#else
typedef float gitr_precision;
netCDF::NcType netcdf_precision = netCDF::ncFloat;
#endif

#ifdef __CUDACC__
    typedef curandState rand_type;
#else
    typedef std::mt19937 rand_type;
#endif


int main(int argc, char **argv, char **envp) {
  CLI::App app{ "!" };
  std::string file_name = "input/gitrInput.cfg";
  app.add_option( "-c", file_name, "config filepath" );
  CLI11_PARSE( app, argc, argv );
  std::cout << "file_name read from stdin: " << file_name << std::endl;
  typedef std::chrono::high_resolution_clock gitr_time;
  auto gitr_start_clock = gitr_time::now();
  class libconfig_string_query query( file_name );
  class use use( query );
  int use_3d_geom = use.get< int >( use::use_3d_geom );
  int bfield_interp = use.get< int >( use::bfield_interp );
  int presheath_interp = use.get< int >( use::presheath_interp );
  int surface_potential = use.get< int >( use::surface_potential );
 int sheath_efield = use.get< int >( use::sheath_efield );
 int ionization = use.get< int >( use::ionization );
 printf("ionization %d \n", ionization);
 int backgroundCollisions = use.get< int >( use::backgroundCollisions );
 int surface_model = use.get< int >( use::surface_model );


 // Set default processes per node to 1
 int ppn = 1;
 // Set default input file string
 std::string inputFile = file_name;
 read_comand_line_args(argc,argv,ppn,inputFile);

int world_rank = 0;
int world_size = 1;

std::string input_path = "input/";
libconfig::Config cfg, cfg_geom;
cfg.setAutoConvert(true);
cfg_geom.setAutoConvert(true);

auto gitr_flags = new Flags(cfg);


gitr_precision background_Z = 0.0, background_amu = 0.0;

 auto performActionIfWorldRankZero = [&](auto action) {
     if (world_rank == 0) { action(); }
 };

 performActionIfWorldRankZero([&](){
     std::cout << "Open configuration file " << input_path << inputFile << std::endl;
     importLibConfig(cfg, inputFile);
     std::string geomFile;
     getVariable(cfg, "geometry.fileString", geomFile);
     std::cout << "Open geometry file " << input_path + geomFile << std::endl;
     importLibConfig(cfg_geom, input_path + geomFile);
     std::cout << "Successfully staged input and geometry file " << std::endl;
 });


 performActionIfWorldRankZero([&](){
     getVariable(cfg, "backgroundPlasmaProfiles.Z", background_Z);
     getVariable(cfg, "backgroundPlasmaProfiles.amu", background_amu);
 });


 int nR_Bfield = 1, nY_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;
 std::string bfieldCfg = "backgroundPlasmaProfiles.Bfield.";
 std::string bfieldFile;

 performActionIfWorldRankZero([&](){
     importVectorFieldNs(cfg, input_path, bfield_interp, bfieldCfg, nR_Bfield, nY_Bfield, nZ_Bfield, bfieldFile);
 });

 sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridy(nY_Bfield), bfieldGridz(nZ_Bfield);
 n_Bfield = nR_Bfield * nY_Bfield * nZ_Bfield;
 sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

 performActionIfWorldRankZero([&](){
     importVectorField(cfg, input_path, bfield_interp, bfieldCfg, nR_Bfield, nY_Bfield, nZ_Bfield,
                       bfieldGridr.front(), bfieldGridy.front(), bfieldGridz.front(),
                       br.front(), by.front(), bz.front(), bfieldFile);
   });

 std::cout << "Start of geometry import" << std::endl;
 int nLines = 1, nSurfaces = 0;
 performActionIfWorldRankZero([&](){
     try {
         libconfig::Setting &geom = cfg_geom.lookup("geom");
         nLines = geom["x1"].getLength();
         std::cout << "nLines " << nLines << std::endl;
     } catch (const libconfig::SettingNotFoundException &nfex) {
         std::cerr << "No 'geom' setting in configuration file." << std::endl;
     }
 });

 sim::Array<Boundary> boundaries(nLines + 1, Boundary());
 performActionIfWorldRankZero([&](){
   nSurfaces = importGeometry(cfg_geom, boundaries, use_3d_geom );
   std::cout << "Starting Boundary Init... nSurfaces " << nSurfaces << std::endl;
 });

  // Background Plasma Temperature Initialization
  int nR_Temp = 1, nY_Temp = 1, nZ_Temp = 1, n_Temp = 1;
  int nR_Dens = 1, nY_Dens = 1, nZ_Dens = 1, n_Dens = 1;
  int nR_flowV = 1, nY_flowV = 1, nZ_flowV = 1, n_flowV = 1;
  n_Temp = nR_Temp * nY_Temp * nZ_Temp;

  sim::Array<gitr_precision> TempGridr(nR_Temp), TempGridz(nZ_Temp), TempGridy(nY_Temp);
  sim::Array<gitr_precision> DensGridr(nR_Dens), DensGridz(nZ_Dens), DensGridy(nY_Dens);

  n_Dens = nR_Dens * nY_Dens * nZ_Dens;
  n_flowV = nR_flowV * nY_flowV * nZ_flowV;

  std::vector<double> flowVr(n_flowV), flowVz(n_flowV), flowVt(n_flowV);

  sim::Array<gitr_precision> ni(n_Dens), ne(n_Dens);
  sim::Array<gitr_precision> ti(n_Temp), te(n_Temp);

// read in the plasma background profiles
std::tuple<int, int, int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> plasmaResult = processPlasmaBackground();
int nTemperatures = std::get<0>(plasmaResult);
int nDensities = std::get<1>(plasmaResult);
int nVelocities = std::get<2>(plasmaResult);

te = std::get<3>(plasmaResult);
ti = std::get<4>(plasmaResult);
ne = std::get<5>(plasmaResult);
ni = std::get<6>(plasmaResult);

flowVr = std::get<7>(plasmaResult);
flowVt = std::get<8>(plasmaResult);
flowVz = std::get<9>(plasmaResult);


// initialize plasma conditions at the surface boundary
 std::for_each(boundaries.begin(), boundaries.end() - 1,
               boundary_init(background_Z, background_amu, ni.data(), ne.data(), br.data(),
                             bz.data(), by.data(), ti.data(), te.data(), surface_potential, use_3d_geom ));

 gitr_precision dt;
 int nP = 0 , nT = 0;
 gitr_precision max_dt = 1.0e5;
 if (world_rank == 0) {
   if (cfg.lookupValue("timeStep.dt", dt) &&
       cfg.lookupValue("timeStep.nT", nT)) {
         std::cout << "Number of time steps: " << nT << " With dt = " << dt << std::endl;
   }
   else
   {
     std::cout << "ERROR: could not get nT, dt, or nP from input file" << std::endl;
   }
 }

   // read in particle source
   std::string ncParticleSourceFile;
   getVariable(cfg, "particleSource.ncFileString", ncParticleSourceFile);
   std::string particleSourceFile = ncParticleSourceFile;
   ParticleData particleData = readParticleData(particleSourceFile);
   nP = particleData.materialName.size();
   long nParticles = nP;
   gitr_precision x, y, z, vx, vy, vz, amu, Z, charge;
   std::string materialName;
   sim::Array<gitr_precision>  px(nP), py(nP), pz(nP), pvx(nP), pvy(nP), pvz(nP),  pZ(nP), pamu(nP), pcharge(nP);
   sim::Array<std::string> pmaterialName(nP);
   auto particleArray = new Particles(nP,1,cfg,gitr_flags);
   initializeParticleArray(particleData, particleArray, px, py, pz, pvx, pvy, pvz, pZ, pamu, pcharge, dt, pmaterialName);

// 
  sim::Array<int> nPPerRank(world_size, 0), pStartIndx(world_size, 0), pDisplacement(world_size, 0), pHistPerNode(world_size, 0), nActiveParticlesOnRank(world_size, 0);
  int countP = 0;
 if (nP >= world_size) {
   for (int i = 0; i < world_size; i++) {
     nPPerRank[i] = std::floor(nP / world_size);
     if (i == 0) {
       nPPerRank[i] = nPPerRank[i] + nP % world_size;
     }
     pStartIndx[i] = countP;
     countP = countP + nPPerRank[i];
     std::cout << "countP " << countP << std::endl;
   }
 } else {
   for (int i = 0; i < nP; i++) {
     nPPerRank[i] = 1;
     pStartIndx[i] = countP;
     countP = countP + nPPerRank[i];
   }
 }

 for (int i = 0; i < world_size; i++) {
   nActiveParticlesOnRank[i] = nPPerRank[i];
 }

 std::uniform_real_distribution<gitr_precision> dist(0, 1e6);
 thrust::counting_iterator<std::size_t> particleBegin(pStartIndx[world_rank]);
 thrust::counting_iterator<std::size_t> particleEnd(pStartIndx[world_rank] + nActiveParticlesOnRank[world_rank] );
 thrust::counting_iterator<std::size_t> particleOne(1);
 thrust::counting_iterator<std::size_t> particleZero(0);
 auto randInitStart_clock = gitr_time::now();


#if USE_CUDA
 sim::Array<rand_type> state1(nParticles);
#else
 sim::Array<rand_type> state1(nParticles);
#endif


 if( ionization > 0 || backgroundCollisions > 0 )
   {
     #if USE_CUDA
       int *dev_i;
       cudaMallocManaged(&dev_i, sizeof(int));
       dev_i[0] = 0;
       thrust::for_each(thrust::device, particleBegin, particleEnd, curandInitialize<>( &state1.front(), 12832 ));
     #else
       std::random_device randDeviceInit;
       thrust::for_each( thrust::device, particleBegin, particleEnd, curandInitialize<>( &state1.front(), 12832 ) );
     #endif
     #if USE_CUDA
       cudaDeviceSynchronize();
     #endif
   }

  pusher pusher0( particleArray, boundaries.data(), nLines, nR_Bfield, nZ_Bfield, bfieldGridr.data(), &bfieldGridz.front(), &br.front(), &bz.front(), &by.front(), gitr_flags);

  // geometryCheck geometryCheck0(particleArray, nLines, &boundaries[0], gitr_flags);

  // backgroundCollisionsMCC backgroundCollisionsMCC0(particleArray, &state1.front(), flowVr[0], flowVz[0], flowVt[0], ne[0], ni[0], ti[0], te[0], background_Z, background_amu);

  // elementaryProcesses<rand_type> elementaryProcesses0( particleArray, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
  //      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
  //      &TempGridz.front(), &te.front(), gitr_flags );

//  surfaceReactions surfaceReactions0(  particleArray, &state1.front(), nLines, &boundaries[0], surfaces, gitr_flags);

 auto start_clock = gitr_time::now();
 std::cout << "Starting main loop" << std::endl;
 #if __CUDACC__
 cudaDeviceSynchronize();
 #endif

  // Diagnostic
     int subSampleFac = 1;
     std::string outputFileName;
     if (world_rank == 0) {
       if (cfg.lookupValue("diagnostics.trackSubSampleFactor", subSampleFac)) {
         std::cout << "Tracks subsample factor imported " << subSampleFac << std::endl;
               // import output file name
         cfg.lookupValue("diagnostics.dumpTrajectories", outputFileName);
       } else {
         std::cout
             << "ERROR: Could not get tracks sub sample info from input file "
             << std::endl;
       }
     }
     // Define the domain
     Domain domain;
     std::string domainBounds = "domain.";
     getVariable(cfg, domainBounds + "xmin", domain.xmin);
     getVariable(cfg, domainBounds + "xmax", domain.xmax);
     getVariable(cfg, domainBounds + "ymax", domain.ymax);
     getVariable(cfg, domainBounds + "ymin", domain.ymin);
     getVariable(cfg, domainBounds + "zmax", domain.zmax);
     getVariable(cfg, domainBounds + "zmin", domain.zmin);

// Main loop
 for (int tt = 0; tt < nT; tt++) {
     thrust::for_each(thrust::device, particleBegin, particleEnd,
         [&] __device__ (auto& particle) {
             pusher0(particle);
            // geometryCheck0(particle);
            //  if (ionization > 0 ) elementaryProcesses0(particle);
            //  if (backgroundCollisions > 0) backgroundCollisionsMCC0(particle);
            //  if (surface_model > 0) surfaceReactions0(particle);

         }

     );
   // Store particle data
   if (tt % subSampleFac == 0)
       storeParticleData(outputFileName, particleArray, nP, domain, tt);
 }

 #if __CUDACC__
 cudaDeviceSynchronize();
 #endif

 auto finish_clock = gitr_time::now();
 std::chrono::duration<gitr_precision> fs = finish_clock - start_clock;
 printf("Time taken is %6.3f (secs) \n", fs.count());
 printf("Time taken per step is %6.3f (secs) \n", fs.count() / (gitr_precision)nT);

return 0;
 
}
