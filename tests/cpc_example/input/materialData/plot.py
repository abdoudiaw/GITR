import netCDF4 as nc

data = nc.Dataset("surfaceReactions_O_on_W.nc", "r")

print(data['rfyld'][:])
