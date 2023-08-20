#import netCDF4 as nc
#
#
#data = nc.Dataset("plasmaProfile.nc","r")
## // get variables and data
#
##print(data.variables.keys())  # dict_keys(['ne', 'te', 'ni', 'ti', 'vx', 'vy', 'vz'])
###print(data.variables['x'][:])
###print(data.variables['y'][:])
###print(data.variables['z'][:])
###print(data.variables['materialName'][:])
###
###data.close()
###


import numpy as np
import netCDF4 as nc

# Open the netCDF file for reading and writing
data = nc.Dataset("plasmaProfile.nc", "r+")

# Check the variables present in the dataset
print(data.variables.keys())  # e.g., dict_keys(['ne', 'te', 'ni', 'ti', 'vx', 'vy', 'vz'])

# Replace the data for the variables
#data.variables['vx'][:] = np.zeros_like(data.variables['vx'][:])
#data.variables['vy'][:] = np.zeros_like(data.variables['vy'][:])
data.variables['vx'][:] = np.full_like(data.variables['vx'][:], -500)
data.variables['vy'][:] = np.full_like(data.variables['vy'][:], -100)
data.variables['vz'][:] = np.full_like(data.variables['vz'][:], -2000)

# Close the dataset after making changes
data.close()
