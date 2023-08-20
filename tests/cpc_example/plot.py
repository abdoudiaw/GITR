import netCDF4 as nc
import numpy  as np
import libconf, io


B= 1
charge = 1.60217662e-19
mass = 1.6726219e-27
amu_oxygen = 16

wc = charge * B / (mass * amu_oxygen)
q_over_m = charge / mass

T = 2*np.pi/ wc

dt = 0.05 * T
print(dt)


#exit()


data = nc.Dataset("O.nc","r")
for var in data.variables:
    print(var)
print(len(data['x'][:]))

##exit()
### // get variables and data
#
##print(data.variables.keys())
#x=data.variables['vx'][:]
#y=data.variables['vy'][:]
#z=data.variables['vz'][:]
#
#
#print(x,y,z)

exit()
#data = nc.Dataset("plasmaProfile.nc", "r")
#for key in data.variables:
#    print(key, data[key][:])
##print(data.variables)
#exit()
filename="input/gitrGeometryPointPlane3d.cfg"
with io.open(filename) as f:
    config = libconf.load(f)

x1 = np.array(config.geom.x1)
x2 = np.array(config.geom.x2)
x3 = np.array(config.geom.x3)
z1 = np.array(config.geom.z1)
z2 = np.array(config.geom.z2)
z3 = np.array(config.geom.z3)
Z = np.array(config.geom.Z)
y1 = np.array(config.geom.y1)
y2 = np.array(config.geom.y2)
y3 = np.array(config.geom.y3)

#import netCDF4 as nc
#import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#
#def check_geometry_closure(x, y, z):
#    vertex_dict = {}
#    for tri_vertices in zip(x, y, z):
#        for vertex in tri_vertices:
#            if vertex not in vertex_dict:
#                vertex_dict[vertex] = []
#        for i in range(3):
#            vertex1 = tri_vertices[i]
#            vertex2 = tri_vertices[(i + 1) % 3]
#            vertex_dict[vertex1].append(vertex2)
#    for vertex, connections in vertex_dict.items():
#        if len(connections) != 2:
#            return False
#    return True
#
#def find_open_edges(x1, x2, x3, y1, y2, y3, z1, z2, z3):
#    tris = np.vstack((x1, x2, x3, y1, y2, y3, z1, z2, z3)).T
#    edges = set()
#    open_edges = set()
#
#    for tri in tris:
#        for i in range(3):
#            edge = (tri[i], tri[(i + 1) % 3])
#            reverse_edge = (tri[(i + 1) % 3], tri[i])
#            if edge in edges:
#                edges.remove(edge)
#            elif reverse_edge in edges:
#                edges.remove(reverse_edge)
#            else:
#                edges.add(edge)
#
#    for edge in edges:
#        open_edges.add(edge[0])
#        open_edges.add(edge[1])
#
#    return open_edges
#
#
#def close_geometry(x1, x2, x3, y1, y2, y3, z1, z2, z3):
#    is_closed = check_geometry_closure(np.concatenate([x1, x2, x3]), np.concatenate([y1, y2, y3]), np.concatenate([z1, z2, z3]))
#
#    if is_closed:
#        print("The geometry is already closed.")
#        return x1, x2, x3, y1, y2, y3, z1, z2, z3
#    else:
#        print("Geometry is open")
#
#
#        open_edges = find_open_edges(x1, x2, x3, y1, y2, y3, z1, z2, z3)
#
#        if len(open_edges) == 0:
#            print("No open edges found in the geometry.")
#            return x1, x2, x3, y1, y2, y3, z1, z2, z3
#
#        new_tris = []
#        for edge in open_edges:
#            v1, v2 = edge
#            new_tri = [v1, v2, v1]
#            new_tris.append(new_tri)
#
#        new_tris = np.array(new_tris).T
#
#        x1 = np.concatenate([x1, new_tris[0]])
#        x2 = np.concatenate([x2, new_tris[1]])
#        x3 = np.concatenate([x3, new_tris[2]])
#        y1 = np.concatenate([y1, new_tris[3]])
#        y2 = np.concatenate([y2, new_tris[4]])
#        y3 = np.concatenate([y3, new_tris[5]])
#        z1 = np.concatenate([z1, new_tris[6]])
#        z2 = np.concatenate([z2, new_tris[7]])
#        z3 = np.concatenate([z3, new_tris[8]])
#
#        return x1, x2, x3, y1, y2, y3, z1, z2, z3
#
#
## Modify the tris to close the geometry
#x1, x2, x3, y1, y2, y3, z1, z2, z3 = close_geometry(x1, x2, x3, y1, y2, y3, z1, z2, z3)


import libconf

filename = "input/gitrGeometryPointPlane3d.cfg"

# Load the configuration from the file
with open(filename) as f:
    config = libconf.load(f)

# Update the geometry values with the corrected ones
config.geom.x1 = x1.tolist()
config.geom.x2 = x2.tolist()
config.geom.x3 = x3.tolist()
config.geom.y1 = y1.tolist()
config.geom.y2 = y2.tolist()
config.geom.y3 = y3.tolist()
config.geom.z1 = z1.tolist()
config.geom.z2 = z2.tolist()
config.geom.z3 = z3.tolist()

# Write the updated configuration back to the file
with open(filename, "w") as f:
    libconf.dump(config, f)


# Plot the tris and particles in 3D as a mesh
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot tris as a mesh
ax.plot_trisurf(x1, y1, z1, linewidth=0.2, antialiased=True, color='blue', alpha=0.5)
ax.plot_trisurf(x2, y2, z2, linewidth=0.2, antialiased=True, color='red', alpha=0.5)
ax.plot_trisurf(x3, y3, z3, linewidth=0.2, antialiased=True, color='green', alpha=0.5)



# Plot particles
ax.scatter(x, y, z, c='m', marker='o')

# Display the plot
plt.show()
