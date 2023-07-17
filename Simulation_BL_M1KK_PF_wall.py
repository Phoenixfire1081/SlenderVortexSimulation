import numpy as np
import sys
import time

sys.path.append('src/')

from RunSimulation import initializeAndRun
from plotTools import plotFunc

start_time = time.time()

#---------------------------------------------------------------------#

# ---------Analysis of Multiscale Data from the Geosciences---------- #

# Author: Abhishek Harikrishnan
# Email: abhishek.harikrishnan@fu-berlin.de
# Last updated: 11-07-2023
# Boundary layer flow test template with M1 KK method and partial FORTRAN mode

#---------------------------------------------------------------------#

# SET ALL PARAMETERS HERE

# ---------------------
# Simulation parameters
# ---------------------

# Number of points
# This should be an odd number. 
NP = 701 

# Time step
ts = 0.05

# Number of time steps
nsteps = 500

# Core size, dttm
epsilon = 0.05 

# Set overlap requirement. If overlapRequirement = epsilon, this 
# satisfies the minimum overlap condition. 
overlapRequirement = epsilon/3

# Enable Local Induction Approximation?
# If _LIA = False, then the simulation is carried out with M1 KK scheme.
_LIA = False

# Set the numerical core coefficient for the M1 KK scheme
# c_ttm = -0.4202 was obtained from Knio 2000
# This is dependent on the core smoothing function and can be evaluated analytically.
c_ttm = -0.4202

# Set additional M1 KK parameters. This is for the method 1 optimization 
# technique of  Knio 2000. K_M1KK and Phi_M1KK should be > 1. 
K_M1KK = 3
Phi_M1KK = 2

# Set circulation value
gamma_param = 0.00001

# Continue from previous calculation?
# If _restart = True, then data is read from the restart files and calculation continues.
_restart = False

# Filename for output
filename = 'R1.txt'

# Write out data every N steps.
writeEveryNSteps = 50

# Enable FORTRAN for better performance?
# NOTE: Works only for the M1 KK method.
# _partialFORTRAN runs only the uttm computation. Time integration is still handled by python
# _fullFORTRAN runs the time integration as well. This is faster than _partialFORTRAN
_partialFORTRAN = True
_fullFORTRAN = False

# Include wall boundary condition?
_image = True

# ---------------------------------
# IC related parameters for hairpin
# ---------------------------------

# Amplitude
A = 10/1000

# Inclination
alpha = 0.1*np.pi/180

# Spread parameter
beta = 7.5*1000

# Starting x distance
xinit = 1/1000

# Starting y distance
yinit = 30/1000

# Spanwise length
L = 200/1000

# Number of periodic boxes for infinitely long filaments
n_b = 8 

# Uses matplotlib to visualize data.
_visualize = True

# --------------------------
# Background flow parameters
# --------------------------

# Uniform flow. Set True/False.
_uniform = False

# If _uniform = True, set velocity.
Ub = 4

# Shear flow. Set True/False.
_shear = False

# If _shear = True, set velocity and height above which flow becomes uniform again.
# If logspaced = True, the nodes are spaced out unevenly with the outer nodes
# spaced far apart than the inner nodes. 
maxV = 30.0
yShear = 1.5
logspaced = False

# Boundary layer flow from data. Set True/False.
_BL = True

# If _BL = True, set friction velocity u_tau, kinematic viscosity nu to calculate the viscous length scale delta_nu.
u_tau = 0.0472
nu = 2/1000000
delta_nu = nu/u_tau

# Additionally, set location of the velocity profiles and the wall-normal height profile. 
if _BL:
	uVelBL = np.loadtxt('../Profiles/FR.02/VelX.txt')/(u_tau*1000)
	vVelBL = np.loadtxt('../Profiles/FR.02/VelY.txt')/(u_tau*1000)
	wVelBL = np.loadtxt('../Profiles/FR.02/VelZ.txt')/(u_tau*1000)
	yPlus = np.loadtxt('../Profiles/y_dir_FR02.txt')/1000
	yPlus = yPlus[:len(wVelBL)]
else:
	uVelBL = np.array([0.0, 0.0, 0.0])
	vVelBL = np.array([0.0, 0.0, 0.0])
	wVelBL = np.array([0.0, 0.0, 0.0])
	yPlus = np.array([0.0, 0.0, 0.0])

# NOTE: if _uniform = False, _shear = False and _BL = False, then there is no (or stagnant) background flow.

# -------------------
# Default parameters.
# -------------------

# Number of vortices. Currently admits only one.
NF = 1 

# Numpy data type. Use np.float64 for best accuracy. np.float32 is faster.
_dtype = np.float64

# Viscosity. This not implemented yet. 
nu_bar_param = 0 

# Axial velocity. This not implemented yet. 
m_0_param = 0 

# Initial core radius. This is used internally.
delta_0_bar_param = epsilon 

# --------------
# Run simulation
# --------------

Ux, Uy, Uz = initializeAndRun(NP, NF, _dtype, ts, u_tau, nu, delta_nu, nsteps, epsilon, nu_bar_param, \
A, alpha, beta, xinit, yinit, L, _uniform, Ub, _shear, maxV, yShear, _BL, \
uVelBL, vVelBL, wVelBL, yPlus, gamma_param, m_0_param, delta_0_bar_param, \
K_M1KK, Phi_M1KK, c_ttm, n_b, _restart, filename, _LIA, overlapRequirement, writeEveryNSteps, \
_partialFORTRAN, _fullFORTRAN, _image, logspaced)

if _partialFORTRAN:
	import os
	os.system('rm -rf fio.dat')

if _fullFORTRAN:
	import os
	os.system('rm -rf fio.dat')
	os.system('mv src/full_fortran/fout.dat .')
	os.system('mv fout.dat ' + filename)

print('Total time taken: ' + str(time.time() - start_time) + ' s')

if _visualize:
	
	plotFunc(filename, nsteps, NP, writeEveryNSteps, 1000)


