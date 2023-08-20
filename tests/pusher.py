def boris_push(x_i,y_i,z_i,  ux_i, uy_i, uz_i, Ex, Ey, Ez, Bx, By, Bz, q_over_m, dt):
    """
    Boris pusher for nonrelativistic plasma
    """

    # Half acceleration by electric field
    ux_minus = ux_i + 0.5 * q_over_m * Ex * dt
    uy_minus = uy_i + 0.5 * q_over_m * Ey * dt
    uz_minus = uz_i + 0.5 * q_over_m * Ez * dt

    # Rotation by magnetic field
    t_x = q_over_m * Bx * 0.5 * dt
    t_y = q_over_m * By * 0.5 * dt
    t_z = q_over_m * Bz * 0.5 * dt

    s_x = 2*t_x / (1 + t_x**2 + t_y**2 + t_z**2)
    s_y = 2*t_y / (1 + t_x**2 + t_y**2 + t_z**2)
    s_z = 2*t_z / (1 + t_x**2 + t_y**2 + t_z**2)

    ux_prime = ux_minus + uy_minus * t_z - uz_minus * t_y
    uy_prime = uy_minus + uz_minus * t_x - ux_minus * t_z
    uz_prime = uz_minus + ux_minus * t_y - uy_minus * t_x

    ux_plus = ux_minus + uy_prime * s_z - uz_prime * s_y
    uy_plus = uy_minus + uz_prime * s_x - ux_prime * s_z
    uz_plus = uz_minus + ux_prime * s_y - uy_prime * s_x

    # Half acceleration by electric field again
    ux_f = ux_plus + 0.5 * q_over_m * Ex * dt
    uy_f = uy_plus + 0.5 * q_over_m * Ey * dt
    uz_f = uz_plus + 0.5 * q_over_m * Ez * dt

    # updaticles positions

    x_f = x_i + ux_f * dt
    y_f = y_i + uy_f * dt
    z_f = z_i + uz_f * dt

    return x_f, y_f, z_f, ux_f, uy_f, uz_f

import numpy as np  
# Path: tests/pusher.py
q = 1
m = 1
dt = 0.05
n_time_steps = 100
n_particles = 1
# initial positions
x = np.zeros(n_particles)
y = np.zeros(n_particles)
z = np.zeros(n_particles)

# initial velocities
ux = np.ones(n_particles)
uy = np.ones(n_particles)
uz = np.zeros(n_particles)

# initial electric field
Ex = np.zeros(n_particles)
Ey = np.zeros(n_particles)
Ez = np.zeros(n_particles)

# initial magnetic field
Bx = np.ones(n_particles)
By = np.zeros(n_particles)
Bz = np.zeros(n_particles)

# run the pusher
positions = np.zeros((n_time_steps, n_particles, 3))
velocities = np.zeros((n_time_steps, n_particles, 3))
for i in range(n_time_steps):
    x, y, z, ux, uy, uz = boris_push(x, y, z, ux, uy, uz, Ex, Ey, Ez, Bx, By, Bz, q/m, dt)
    positions[i, :, 0] = x
    positions[i, :, 1] = y
    positions[i, :, 2] = z

    velocities[i, :, 0] = ux
    velocities[i, :, 1] = uy
    velocities[i, :, 2] = uz


    # print(x,y,z,ux,uy,uz)

positions = np.reshape(positions, (n_time_steps, n_particles*3))
positions=np.array(positions)
print(positions)

velocities = np.reshape(velocities, (n_time_steps, n_particles*3))
velocities=np.array(velocities)
print(velocities)
