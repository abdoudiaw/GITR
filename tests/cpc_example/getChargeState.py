# import matplotlib.pyplot as plt
# from collections import Counter

# # Function to parse dump.step file and extract y and z values
# def extract_data(filename):
#     z_values = []
#     with open(filename, 'r') as f:
#         lines = f.readlines()
        
#         # Initializing flag for data extraction
#         extract_data = False
        
#         for line in lines:
#             # If the line indicates atomic data, set the flag to start data extraction
#             if 'ATOMS id' in line:
#                 extract_data = True
#                 continue
#             if extract_data:
#                 data = line.strip().split()
#                 # Extract y and z values
#                 # skip first line 
#                 # get the values from the 9th column
#                 z_values.append(float(data[8]))
# #                print("z_values", z_values)
                
#                 # Reset the flag as we have read the required data for that timestep
#                 extract_data = True
#     return z_values


# # Read data
# # Y = []
# charge_data = []
# r_gyro =  0.001757002937925605
# for step in range(0,   99900 ,100):
#     z = extract_data("output/dump."+str(step))
#     # print(z)
#     charge_data.append(z)
# # Initialize a dictionary to keep track of the charge count over time for each species.
# charge_count_over_time = {}

# # Loop over each time step to count the occurrences of each charge state.
# for t, charges_at_t in enumerate(charge_data):
#     charge_count = Counter(charges_at_t)
    
#     for charge, count in charge_count.items():
#         if charge not in charge_count_over_time:
#             charge_count_over_time[charge] = []
        
#         # Extend the count list to reach the current time step if needed.
#         while len(charge_count_over_time[charge]) <= t:
#             charge_count_over_time[charge].append(0)
        
#         charge_count_over_time[charge][t] = count

# dt = 1 
# # Plotting
# for charge, counts in charge_count_over_time.items():
#     plt.plot(range(len(counts)), counts, label=f"Charge {charge}")

# plt.xlabel("Time step [ns]")
# plt.ylabel("Charge state")
# plt.title("Number of each charge state over time")
# plt.legend()
# plt.show()


import numpy as np
import matplotlib.pyplot as plt

# Parameters
time_steps = 10000  # number of time steps
dt = 0.1  # time step size
Z_max = 8  # maximum Z value
n_e = 1.0  # electron density

# Initialize arrays for n_Z, S_Z, and alpha_Z
n_Z = np.zeros(Z_max + 1)  # n_Z for Z=0 to Z=Z_max
S_Z = np.ones(Z_max + 1) * 0.1  # rate constants S_Z, here assumed constant for simplicity
alpha_Z = np.ones(Z_max + 1) * 0.05  # rate constants alpha_Z, here assumed constant for simplicity

# Initial condition: Assume all ions are in Z=0 state
n_Z[0] = 1.0

# Time evolution
time_series = np.zeros((time_steps, Z_max + 1))
for t in range(time_steps):
    time_series[t, :] = n_Z  # Store current state
    
    # Calculate rate of change (dn_Z/dt) using Euler's method
    dn_Z_dt = np.zeros(Z_max + 1)
    
    # Boundary conditions: Z=0 and Z=Z_max
    dn_Z_dt[0] = n_e * (-n_Z[0] * S_Z[0] + n_Z[1] * alpha_Z[1])
    dn_Z_dt[Z_max] = n_e * (n_Z[Z_max - 1] * S_Z[Z_max - 1] - n_Z[Z_max] * alpha_Z[Z_max])
    
    # Loop over other Z values
    for Z in range(1, Z_max):
        dn_Z_dt[Z] = n_e * (n_Z[Z - 1] * S_Z[Z - 1] - n_Z[Z] * S_Z[Z] + n_Z[Z + 1] * alpha_Z[Z + 1] - n_Z[Z] * alpha_Z[Z])
    
    # Update n_Z
    n_Z += dn_Z_dt * dt

# Plotting
plt.figure(figsize=(10, 6))
for Z in range(Z_max + 1):
    plt.plot(np.arange(time_steps) * dt, time_series[:, Z], label=f"Z={Z}")

plt.xlabel("Time")
plt.ylabel("n_Z")
plt.legend()
plt.show()
