import netCDF4 as nc


data = nc.Dataset('particleSource.nc')

print(data.variables.keys())
charge = data.variables['charge'][:]

print(charge)

IonizationState = data.variables['IonizationState'][:]

print(IonizationState)
