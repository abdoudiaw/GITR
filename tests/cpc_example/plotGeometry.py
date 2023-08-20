import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits import mplot3d
import io, libconf
import os

#######################--------------------#######################
#######################--------------------#######################

filename="input/gitrGeometryPointPlane3d.cfg"
#filename="gitrGeometry.cfg"

with io.open(filename) as f:
    config = libconf.load(f)
    
data=config['geom']
coords=['x1','x2','x3','y1','y2','y3','z1','z2','z3']
[x1,x2,x3,y1,y2,y3,z1,z2,z3]=[data.get(var) for var in coords]
 
 
 
print('x', min(x1)*1., min(x2)*1.e-0)
print('y', min(y1)*1.e-0, min(y2)*1.e-0)
print('z', min(z1)*1.e-0, min(z2)*1.e-0)

print('x', max(x1)*1., max(x2)*1.e-0)
print('y', max(y1)*1.e-0, max(y2)*1.e-0)
print('z', max(z1)*1.e-0, max(z2)*1.e-0)
#exit(0)
import numpy as np

# Defining the points
A = np.array([0.01501, 0.00856])
B = np.array([-0.015, 0.00856])
C = np.array([-0.015, -0.007])


print(-np.sin(2*np.pi/180), np.cos(2*np.pi/180))

#exit()
# Calculate the lengths of the sides
a = np.linalg.norm(B - C)  # BC
b = np.linalg.norm(A - C)  # AC
c = np.linalg.norm(A - B)  # AB

# Use the law of cosines to calculate the angles
cos_angle_A = (b**2 + c**2 - a**2) / (2 * b * c)
cos_angle_B = (a**2 + c**2 - b**2) / (2 * a * c)
cos_angle_C = (a**2 + b**2 - c**2) / (2 * a * b)

angle_A = np.arccos(cos_angle_A) * (180 / np.pi)  # Convert to degrees
angle_B = np.arccos(cos_angle_B) * (180 / np.pi)
angle_C = np.arccos(cos_angle_C) * (180 / np.pi)

print(f"Angle A: {angle_A:.2f}°")
print(f"Angle B: {angle_B:.2f}°")
print(f"Angle C: {angle_C:.2f}°")

# Plotting the points A, B, and C
plt.scatter(*A, color='red', label='A')
plt.scatter(*B, color='blue', label='B')
plt.scatter(*C, color='green', label='C')

# Labeling the points
plt.text(A[0], A[1], 'A', fontsize=12, ha='right')
plt.text(B[0], B[1], 'B', fontsize=12, ha='right')
plt.text(C[0], C[1], 'C', fontsize=12, ha='right')
plt.plot(x3,z3, 'g^')

plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.legend(loc="upper left")
#plt.show()
plt.show()
#exit()
#import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#
## Assuming the previous code with x1, x2, x3, y1, y2, y3, z1, z2, z3 has been run
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#colors = plt.cm.viridis(np.linspace(0, 1, len(x1)))  # Generate distinct colors
#
#for i in range(len(x1)):
#    # Extracting the three points for the current plane
#    p1 = [x1[i], y1[i], z1[i]]
#    p2 = [x2[i], y2[i], z2[i]]
#    p3 = [x3[i], y3[i], z3[i]]
#
#    # Plotting each triangle (plane) with a separate color
#    verts = [np.array([p1, p2, p3])]
#    ax.add_collection3d(Poly3DCollection(verts, color=colors[i], alpha=0.5))
#
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')
#
#plt.show()
#
#exit()
 
#import numpy as np
#
## Extracting points
#p1 = np.array([x1[0], y1[0], z1[0]])
#p2 = np.array([x2[0], y2[0], z2[0]])
#p3 = np.array([x3[0], y3[0], z3[0]])
#
## Compute vectors a and b
#a = p2 #- p1
#b = p3 #- p1
#
## Compute the normal to the plane
#n = np.cross(a, b)
#
## Normalize the normal
#n = n / np.linalg.norm(n)
#
## Vector v directed towards z
#v = np.array([1, 0, 0])
#
## Compute the dot product
#dot = np.dot(n, v)
#
## Compute the angle
#theta = np.arccos(dot)
#
## Convert the angle from radians to degrees
#angle_in_degrees = np.degrees(theta)
#print("The angle between the vector v and the plane is:", angle_in_degrees, "degrees")
#

#exit()
##
##print(x1,x2,x3,y1,y2,y3,z1,z2,z3)
##exit()
figure = plt.figure(figsize=(14, 10))
ax = figure.add_subplot(111, projection='3d')
ax.autoscale(tight=False)

ax.xaxis.set_tick_params(labelsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.zaxis.set_tick_params(labelsize=20)

# Creating plot
#ax.scatter3D(x, y, z, color = "green")
#ax.plot(np.append(x1,x1[0]), np.append(y1,y1[0]), np.append(z1,z1[0]) ,'bs')
#ax.plot(np.append(x2,x2[0]), np.append(y2,y2[0]), np.append(z2,z2[0]), 'ro')
#ax.plot(np.append(x3,x3[0]), np.append(y3,y3[0]), np.append(z3,z3[0]), 'g^')
ax.plot(x3,y3,z3, 'g^')

# fig, ax = plt.subplots(figsize = (12, 7))

#Z_=data.get('Z')
#for i
#x_direct = Z_/np.asarray(74.0)*0.001
#s=where(Z)
#ax.quiver(x1, y1, z1, x_direct, x_direct,x_direct,linewidths = (1,), edgecolor="cyan");
 # ax.set_title('Quiver plot with one arrow')


#ax.set_title('Final W Particle Positions')
ax.set_zlabel('Z[m]',fontsize=20, labelpad=30)
ax.set_xlabel('X[m]',fontsize=20, labelpad=20)
ax.set_ylabel('Y[m]',fontsize=20, labelpad=20)
#ax.tick_params(axis='z', pad=22)
#ax.set_zlim(0., 50e-3)
#ax.set_xlim(-25e-3, 25e-3)
#ax.set_ylim(-25e-3, 25e-3)
ax.grid(False)
plt.tight_layout()
plt.savefig('geometry.png')
plt.show()

