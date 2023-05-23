import numpy as np
import matplotlib.pyplot as plt
import os

#---------------------------------------------------------------------#

# ---------Analysis of Multiscale Data from the Geosciences---------- #

# Author: Abhishek Harikrishnan
# Email: abhishek.harikrishnan@fu-berlin.de
# Last updated: 19-05-2023
# Plots the velocity profiles for the Ekman flow cases used in
# On the motion of hairpin filaments in the atmospheric boundary layer
# Preprint: https://arxiv.org/abs/2303.09302

#---------------------------------------------------------------------#

def readData(*args):
	
	fullPath = os.path.join(*args)
	return np.loadtxt(fullPath)


# Path for profiles
filepath = 'Profiles/'

# List all files
velocity_files = ['VelX.txt', 'VelY.txt', 'VelZ.txt']

# List Ekman flow cases
flow_cases = ['FR.02', 'FR.07', 'FR.inf'] # FR is the Froude number

# Friction velocity (u_{\tau}) for all cases
frictionFR02 = 0.0472
frictionFR07 = 0.0482
frictionFRinf = 0.0528

# Get velocity data and calculate u^+ = u/u_{\tau}
# FR02

FR02Velx = readData(filepath, flow_cases[0], velocity_files[0]) / frictionFR02
FR02Vely = readData(filepath, flow_cases[0], velocity_files[1]) / frictionFR02
FR02Velz = readData(filepath, flow_cases[0], velocity_files[2]) / frictionFR02

# FR07
FR07Velx = readData(filepath, flow_cases[1], velocity_files[0]) / frictionFR07
FR07Vely = readData(filepath, flow_cases[1], velocity_files[1]) / frictionFR07
FR07Velz = readData(filepath, flow_cases[1], velocity_files[2]) / frictionFR07

# FRinf
FRinfVelx = readData(filepath, flow_cases[2], velocity_files[0]) / frictionFRinf
FRinfVely = readData(filepath, flow_cases[2], velocity_files[1]) / frictionFRinf
FRinfVelz = readData(filepath, flow_cases[2], velocity_files[2]) / frictionFRinf

# Get y-coordinate data. Unlike the horizontal directions, this is not evenly spaced.
# Also, unlike velocity data, this is already in viscous units. 
# To convert back to outer units, y = y^+ (\nu) / u_{\tau}
# Note: \nu = 2/1000000
yp_FR02 = readData(filepath, 'y_dir_FR02.txt')
yp_FR07 = readData(filepath, 'y_dir_FR07.txt')
yp_FRinf = readData(filepath, 'y_dir_FRinf.txt')

# Plot the profiles
fig, (ax1, ax2, ax3) = plt.subplots(1, 3,)
ax1.semilogy(FR02Velx, yp_FR02, 'k--', label = 'S_1')
ax1.semilogy(FR07Velx, yp_FR07, 'k:', label = 'S_2')
ax1.semilogy(FRinfVelx, yp_FRinf, 'k-', label = 'N')
ax1.set_xlabel(r'$u^+$')
ax1.set_ylabel(r'$y^+$')
ax1.set_ylim([1, 1500])
ax1.legend()

ax2.semilogy(FR02Vely, yp_FR02, 'k--', label = 'S_1')
ax2.semilogy(FR07Vely, yp_FR07, 'k:', label = 'S_2')
ax2.semilogy(FRinfVely, yp_FRinf, 'k-', label = 'N')
ax2.set_xlabel(r'$v^+$')
ax2.set_ylim([1, 1500])

ax3.semilogy(FR02Velz, yp_FR02, 'k--', label = 'S_1')
ax3.semilogy(FR07Velz, yp_FR07, 'k:', label = 'S_2')
ax3.semilogy(FRinfVelz, yp_FRinf, 'k-', label = 'N')
ax3.set_xlabel(r'$w^+$')
ax3.set_ylim([1, 1500])

plt.show()
