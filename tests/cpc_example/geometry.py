import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits import mplot3d
import io, libconf
import os

#######################--------------------#######################
#######################--------------------#######################

filename="input/gitrGeometry.cfg"

with io.open(filename) as f:
    config = libconf.load(f)
    
data=config['geom']
coords=['x1','x2','x3','y1','y2','y3','z1','z2','z3']
[x1,x2,x3,y1,y2,y3,z1,z2,z3]=[data.get(var) for var in coords]
    
fig = plt.figure(figsize = (10, 8))
ax = fig.add_subplot(111, projection='3d')

ax.autoscale(tight=False)

ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.zaxis.set_tick_params(labelsize=20)

# Creating plot
#ax.scatter3D(x, y, z, color = "green")
ax.plot(np.append(x1,x1[0]), np.append(y1,y1[0]), np.append(z1,z1[0]) ,'bs')
#ax.plot(np.append(x2,x2[0]), np.append(y2,y2[0]), np.append(z2,z2[0]), 'bo')
#boundaryVector[wallHit].x, 0.0183633
#boundaryVector[wallHit].x, 0.0408485
#boundaryVector[wallHit].x, 0.0919991
A=np.array([0.0183633,0.0408485,0.0919991])
# fig, ax = plt.subplots(figsize = (12, 7))
ax.plot(*A,'ro')
print('x', min(x1)*1.e-3, min(x2)*1.e-3)
print('y', min(y1)*1.e-3, min(y2)*1.e-3)
print('z', min(z1)*1.e-3, min(z2)*1.e-3)


plt.show()

#exit()
print("keys", data.keys())

print(" Max x1: ", max(x1) , " Min x1: ", min(x1))
print(" Max x2: ", max(x2) , " Min x2: ", min(x2))
print(" Max x3: ", max(x3) , " Min x3: ", min(x3))
print(" Max y1: ", max(y1) , " Min y1: ", min(y1))
print(" Max y2: ", max(y2) , " Min y2: ", min(y2))
print(" Max y3: ", max(y3) , " Min y3: ", min(y3))
print(" Max z1: ", max(z1) , " Min z1: ", min(z1))
print(" Max z2: ", max(z2) , " Min z2: ", min(z2))
exit()
Z=  data.get('Z')
print("Z: ", Z)
inDir = data.get('inDir')
surface = data.get('surface')
print("inDir: ", inDir)
print("surface: ", surface)
exit()
length_of_x = len(x1)

# # Add 'surface' and 'inDir' keys to the 'geom' section of the config
# data['surface'] = [1] * length_of_x
# data['inDir'] = [-1] * length_of_x
# data['Z'] = [74] * length_of_x

# # Save the updated configuration back to the file
# with io.open(filename, 'w') as f:
#     libconf.dump(config, f)
    
# print("keys", data.keys())
