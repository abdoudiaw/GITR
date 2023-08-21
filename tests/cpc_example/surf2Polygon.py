import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from mpl_toolkits import mplot3d
import io, libconf
import os
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#######################--------------------#######################
#######################--------------------#######################

filename="input/gitrGeometry.cfg"

with io.open(filename) as f:
    config = libconf.load(f)
    
data = config['geom']
coords = ['x1','x2','x3','y1','y2','y3','z1','z2','z3']
[x1, x2, x3, y1, y2, y3, z1, z2, z3] = [data.get(var) for var in coords]

# Assuming the data represents multiple triangles, and each xi, yi, zi is a list of points

# Form the triangles
triangles = []
for i in range(len(x1)):
    triangle = [(x1[i], y1[i], z1[i]), (x2[i], y2[i], z2[i]), (x3[i], y3[i], z3[i])]
    triangles.append(triangle)



px, py, pz=0.0092, 0.04, 0.04
particle_position = np.array([px, py, pz])


def is_point_in_plane(p, a, b, c):
    # Compute the normal of the plane defined by the triangle
    normal = np.cross(b-a, c-a)
    normal /= np.linalg.norm(normal)
    
    # Check if the point p lies in this plane
    return abs(np.dot(normal, p-a)) < 1e-6  # Using a small threshold to account for numerical precision

def point_inside_3D_triangle(p, a, b, c):
    # First, check if point is in the plane of the triangle
    if not is_point_in_plane(p, a, b, c):
        print("not inside")
        return False

    # If yes, check if it's inside the triangle
    return point_in_triangle(p, a, b, c)

is_inside_any_triangle = False

for i in range(len(x1)):
    a = np.array([x1[i], y1[i], z1[i]])
    b = np.array([x2[i], y2[i], z2[i]])
    c = np.array([x3[i], y3[i], z3[i]])

    if point_inside_3D_triangle(particle_position, a, b, c):
        is_inside_any_triangle = True
        print("it is in triangle")
        break  # This break should be indented to be inside the for loop



exit()

# Plotting the triangles
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for triangle in triangles:
    triangle_np = np.array(triangle)
    
    # Check if we have at least 3 unique points
    unique_points = {tuple(row) for row in triangle_np}
    if len(unique_points) >= 3:
        polygon = Poly3DCollection([triangle_np])
        ax.add_collection3d(polygon)

ax.scatter(*particle_position, 'ro')
#ax.set_xlim(0,  0.0919991)
#ax.set_ylim(0,  0.0919991)
#ax.set_zlim(0,  0.0919991)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
