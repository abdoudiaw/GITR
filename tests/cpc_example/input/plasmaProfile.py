import netCDF4 as nc


data = nc.Dataset("plasmaProfile.nc", "r")

for var in data.variables.keys():
        print(var, data[var][:])
        
