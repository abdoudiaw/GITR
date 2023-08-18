#include "particles.h"
#include "array.h"
#include <iostream>
#include <fstream>


// Domain geometry struct
struct Domain {
    gitr_precision xmin, xmax, ymin, ymax, zmin, zmax;
};

// Function to store particle data in a file
void storeParticleData(const std::string& base_filename, Particles *particleArray, int nP, const Domain& domain, int step_count){

    std::string filename = base_filename + std::to_string(step_count);
    // Open the file in append mode
    // std::ofstream f("output/" + filename, std::ios::out);
// print we are here 
    std::ofstream f("output/" + filename, std::ios::app);
    if (!f) {
        std::cout << "Error opening file: " << filename << std::endl;
        return;
    }
    // Get the number of species using unique materialName values
    std::vector<std::string> materialNames;
    for (int i = 0; i < nP; i++) {
        if (std::find(materialNames.begin(), materialNames.end(), particleArray->materialName[i]) == materialNames.end()) {
            materialNames.push_back(particleArray->materialName[i]);
            // std::cout << "Found material: " << particleArray->materialName[i] << std::endl;
        }
    }
    // get unique material names
    int nSpecies = materialNames.size();

    // create index for each species starting from 1 to nSpecies
    std::map<std::string, int> materialIndex;
    for (int i = 0; i < nSpecies; i++) {
        materialIndex[materialNames[i]] = i+1;
    }    
    // Write the header
    f << "ITEM: TIMESTEP\n";
    f << step_count << '\n';
    f << "ITEM: NUMBER OF ATOMS\n";
    f << nP << '\n';
    f << "ITEM: BOX BOUNDS pp pp pp\n";
    // Write max and min of domain
    f << domain.xmin << ' ' << domain.xmax << '\n';
    f << domain.ymin << ' ' << domain.ymax << '\n';
    f << domain.zmin << ' ' << domain.zmax << '\n';
    f << "ITEM: ATOMS id type x y z vx vy vz\n";
    // Write the particle data to the file: id type x y z vx vy vz where type is the material index

    for (int i = 0; i < nP; i++) {
      f << i+1  << "  " << materialIndex[particleArray->materialName[i]] << "  " 
                << particleArray->x[i] << "  " << particleArray->y[i]
                << "  " << particleArray->z[i] <<  "  " << particleArray->vx[i]
                << "  " << particleArray->vy[i] << "  " << particleArray->vz[i] <<  std::endl;
                // << "  " << particleArray->charge[i] << "  " << particleArray->Z[i]
                // << "  " << particleArray->amu[i] <<  " " << particleArray->materialName[i] <<  std::endl;

    }

    f.close();
}
