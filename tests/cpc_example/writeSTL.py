import numpy as np
import libconf
import io

filename = "input/gitrGeometryPointPlane3d.cfg"

# Load config
with io.open(filename) as f:
    config = libconf.load(f)

data = config['geom']
coords = ['x1', 'x2', 'x3', 'y1', 'y2', 'y3', 'z1', 'z2', 'z3']
[x1, x2, x3, y1, y2, y3, z1, z2, z3] = [data.get(var) for var in coords]

# Convert coordinates into triangular vertices
triangles = np.array([
    [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]],
])

# Write triangles to STL
def write_stl(filename, triangles):
    with open(filename, 'w') as f:
        f.write("solid geometry\n")
        for triangle in triangles:
            f.write("  facet normal 0 0 0\n")  # Assuming no normals for simplicity
            f.write("    outer loop\n")
            for vertex in triangle:
                f.write("      vertex {0} {1} {2}\n".format(*vertex))
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write("endsolid geometry\n")

write_stl('output.stl', triangles)

