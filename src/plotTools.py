import matplotlib.pyplot as plt
import numpy as np

#---------------------------------------------------------------------#

# ---------Analysis of Multiscale Data from the Geosciences---------- #

# Author: Abhishek Harikrishnan
# Email: abhishek.harikrishnan@fu-berlin.de
# Last updated: 02-06-2023
# Plotting functions

#---------------------------------------------------------------------#

def plotFunc(filename, nsteps, NP, writeEveryNSteps):
	
	data = np.loadtxt(filename)
	
	ctrB = 3
	ctrE = NP
	
	ax = plt.figure().gca(projection='3d')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	
	nts = int(nsteps/writeEveryNSteps)
	
	for i in range(nts):
	
		if i == 0:
			# First time step is colored blue
			ax.plot(data[ctrB:ctrE, 0], data[ctrB:ctrE, 1], data[ctrB:ctrE, 2], 'b')
		elif i == nts-1:
			# Last time step is colored red
			ax.plot(data[ctrB:ctrE, 0], data[ctrB:ctrE, 1], data[ctrB:ctrE, 2], 'r')
		else:
			# Rest are colored black
			ax.plot(data[ctrB:ctrE, 0], data[ctrB:ctrE, 1], data[ctrB:ctrE, 2], 'k')
	
		ctrB = ctrB + NP;
		ctrE = ctrE + NP;
	
	ax.grid(False)
	ax.azim = -90
	ax.elev = 90
	plt.show()
