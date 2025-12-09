"""
Plot Standard radial type trajectory against Golden angle radial

Choose Nr=100 as example and plot the first 5 spokes

"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# Define Transformation equations

def transform(x_prime, y_prime, theta):
    x = x_prime*np.cos(theta) - y_prime*np.sin(theta)
    y = x_prime*np.sin(theta) + y_prime*np.cos(theta)
    arr = np.vstack((x, y))
    return arr

x_vals = np.linspace(-100, 100, 200)
y_vals = np.linspace(0, 0, 200)
Nr = 100
theta_standard = np.deg2rad(180/Nr)
theta_golden = np.deg2rad(111.25)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(x_vals, y_vals, 'r', linewidth=0.8)
for i in range(1, 6):
    theta = i * theta_standard
    x, y = transform(x_vals, y_vals, theta)
    ax1.plot(x, y, 'r')

ax1.set_aspect('equal')
ax1.set_title('Standard Radial Trajectory\n(Δθ = 1.8°)', fontsize=18)
ax1.set_xlim([-100, 100])
ax1.set_ylim([-100, 100])
ax1.set_xlabel('x', fontsize=16)
ax1.set_ylabel('y', fontsize=16)
#ax1.set_xticks(fontsize=14)
#ax1.set_yticks(fontsize=14)
ax1.grid(True)
ax1.legend()


ax2.plot(x_vals, y_vals,'r',  linewidth=0.8)
for i in range(1, 5):
    theta = i * theta_golden
    x, y = transform(x_vals, y_vals, theta)
    ax2.plot(x, y, 'r')

ax2.set_aspect('equal')
ax2.set_title('Golden-Angle Radial Trajectory\n(Δθ = 111.25°)', fontsize=18)
ax2.set_xlabel('x', fontsize=18)
ax2.set_ylabel('y', fontsize=18)
ax2.set_xlim([-100, 100])
ax2.set_ylim([-100, 100])
#ax2.set_xticks(fontsize=14)
#ax1.set_yticks(fontsize=14)
ax2.grid(True)
ax2.legend()
ax1.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
plt.suptitle('Comparison of First 5 Spokes of Standard vs Golden-Angle Radial Trajectories', fontsize=20)
#plt.tight_layout()
plt.show()