import numpy as np
import matplotlib.pyplot as plt

# Sample data
#data = """
#-0.00363515 0.00750564 -0.00209876 730 60
#-0.00363515 0.00750564 -0.00209876 731 60
#-0.0035778 0.00375877 -0.00206565 732 60
#-0.00173155 0.00563221 -0.000999714 733 60
#-0.0035778 0.00375877 -0.00206565 734 60
#-0.0016742 0.00188533 -0.000966603 735 60
#"""

# Convert the sample data into a NumPy array
#data = np.array([[float(j) for j in i.split()] for i in data.split('\n') if i])
data = np.loadtxt("surface_positions.txt", unpack=True)
# Extract the columns
x = data[0]
y = data[1]
z = data[2]
surface = data[3].astype(int)
angle = data[4].astype(int)

print(np.unique(angle))

exit()
# Define colors based on angle
color_map = {0: 'red', 60: 'green', 90: 'yellow'}

colors = [color_map[a] for a in angle]

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter plot
sc = ax.scatter(x, y, z, c=colors, depthshade=True)

# Add a legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'{key}°',
                          markersize=10, markerfacecolor=color_map[key]) for key in color_map]
ax.legend(handles=legend_elements)

# Display the plot
plt.show()

