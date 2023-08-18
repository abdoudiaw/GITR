#  ADAS Files
#Introduction
PMI expects an ADAS file for each nuclear charge in the system. These files provide atomic and plasma data necessary for PMI simulations to properly calculate ionization and recombination.

# File Naming
The name of each ADAS file should follow the format: ```ruby  ADAS_RATES_{SPECIES NUCLEAR CHARGE}.nc  ```, where  ```{SPECIES NUCLEAR CHARGE} ``` refers to the nuclear charge of the species.

# File Format
The ADAS files are in the NETCDF3 format and follow the NETCDF3_CLASSIC data model.

# Structure of the File
The file contains the following dimensions:
```ruby 
-  n_Temperatures_Ionize 
-  n_Densities_Ionize 
-  n_Temperatures_Recombine 
-  n_Densities_Recombine 
-  n_ChargeStates_Ionize 
-  n_ChargeStates_Recombine 
```
 
The variables in the file include:
```ruby 
- Atomic_Number (int32)
-  gridTemperature_Ionization (float64)
-  gridDensity_Ionization (float64)
-  gridTemperature_Recombination (float64)
-  gridDensity_Recombination (float64)
-  gridChargeState_Ionization (float64)
-  gridChargeState_Recombination (float64)
-  IonizationRateCoeff (float64)
-  RecombinationRateCoeff (float64)
```
# Storage
Store these ADAS files in the input/adasData directory.

