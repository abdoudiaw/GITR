import netCDF4 as nc


data = nc.Dataset("particleSource.nc","r")
# // get variables and data

print(data.variables.keys())
print(data.variables['x'][:])
print(data.variables['y'][:])
print(data.variables['z'][:])
print(data.variables['materialName'][:])

data.close()

