import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_distribution(alpha, B, file_path, save_path=None):
    # Load data
    data = np.loadtxt(file_path, unpack=True, skiprows=0)
    vx = data[3]
    vy = data[4]
    vz = data[5]

    print("Number of particles reaching the surface are", len(vx))
    
    tanTheta = np.sqrt(vz**2 + vy**2) / abs(vx)
    theta = np.arctan(tanTheta) * 180. / np.pi

    incident_energy_hit = np.sqrt(vx**2 + vy**2 + vz**2)
    
    num_bins = 256
    max_energy = max(25, 1.5 * incident_energy_hit.max())
    
    bin_energy = np.linspace(0., max_energy, num_bins // 2)
    bin_angle = np.linspace(0., 90, num_bins // 2)

    # Calculate the histogram
    heights, xedges, yedges = np.histogram2d(theta, incident_energy_hit, bins=[bin_angle, bin_energy], density=False)

    # Normalize
    heights = heights / heights.max()

    # Plotting
    fig, ax = plt.subplots()
    c = ax.pcolormesh(xedges, yedges, heights.T, shading='auto')
    plt.colorbar(c)
    plt.title("B [T] = " +str(B), "$\\alpha [^\\circ] $ = " + str(alpha), fontsize=16)
    plt.xlabel("$\\theta [^\\circ]$", fontsize=16)
    plt.ylabel("$E / T_e$", fontsize=16)

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.show()

# Test for Z = 8
Z = 8.0
alphas = [0.5, 5.0] 
for alpha in alphas:
    file_path = f"../tests/Z_O_8/distributionZmin_alpha_{alpha}.out"
    save_file_path = f"f_E_O_Z_{Z}_alpha_{alpha}.png"
    plot_distribution(alpha, Z, file_path, save_file_path)
