# PMI Surface Reactions Data

This script allows you to create a surface model file for PMI. PMI requires a surface reactions file to describe the interactions between each species in the system and the material, as well as interactions between the wall material and the wall itself.

# File Naming
The name of the surface reactions file should follow the format: ```ruby surfaceReactions_{SPECIES}_on_{SPECIES}.nc```, where ```ruby{SPECIES}``` represents the name of the plasma species. Please use the periodic table notations, such as "O" for oxygen and "W" for tungsten. Avoid using lowercase or full names like "w" or "tungsten".

# File Format
The surface reactions file follows the NETCDF4 data model and is stored in the HDF5 file format.

# Structure of the File
The file contains the following dimensions:

```ruby
nE 
nA 
nEdistBins 
nEdistBinsRef 
nAdistBins 
The variables in the file include:

spyld (float64)
rfyld (float64)
E (float64)
A (float64)
cosXDist (float64)
cosYDist (float64)
cosXDistRef (float64)
cosYDistRef (float64)
energyDist (float64)
energyDistRef (float64)
eDistEgrid (float64)
eDistEgridRef (float64)
phiGrid (float64)
thetaGrid (float64)
```
# Storage
Store these surface reactions files in the input/materialData directory.
