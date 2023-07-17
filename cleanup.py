import os

#---------------------------------------------------------------------#

# ---------Analysis of Multiscale Data from the Geosciences---------- #

# Author: Abhishek Harikrishnan
# Email: abhishek.harikrishnan@fu-berlin.de
# Last updated: 05-07-2023
# Simple script to clean up main directory

#---------------------------------------------------------------------#

ignoreList = ['LICENSE', 'plotProfiles.py', 'README.md', 'Simulation_stagnant_LIA.py', \
'Simulation_stagnant_M1KK.py', 'cleanup.py', 'Plots', 'Profiles', 'src', '.git',\
'Simulation_stagnant_M1KK_PF.py', 'Simulation_stagnant_M1KK_FF.py', 'Simulation_simpleshear_M1KK_FF.py',\
'Simulation_BL_M1KK_FF_nowall.py']
fileList = os.listdir('.')

for i in fileList:
	if not i in ignoreList:
		os.system('rm -rf ' + i)
