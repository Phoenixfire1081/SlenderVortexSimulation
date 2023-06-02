import numpy as np
import os
from numba import jit, prange

# This is a copy of LIAM1functions which acts as a wrapper for the
# full FORTRAN code.

# C_V is the circumferential component of the 
# core constant which is determined from the 
# similar core solutions of Ting Klein 1991. 
# NOTE: C_W is not calculated since axial flux is
# set to 0.
def C_V(delta_bar):
	return (((1 + 0.577 - np.log(2)) / 2) - np.log(delta_bar))

# This function writes the necessary input data for the partial FORTRAN code.
def fortran_writer(NP, delta_ttm, sigma_1, sigma_2, gamma_param, j, ds, n_b, Ux, Uy, \
	Uz, Ux_s, Uy_s, Uz_s, SIGMA, _image, _uniform, Ub, _shear, maxV, yShear, _BL, yPlus,\
	uVelBL, vVelBL, wVelBL, dt):
	
	fw = open('fio.dat', 'a')
	
	fw.write(str(NP) + '\n')
	fw.write(str(delta_ttm) + '\n')
	fw.write(str(sigma_1) + '\n')
	fw.write(str(sigma_2) + '\n')
	fw.write(str(gamma_param[j]) + '\n')
	fw.write(str(ds) + '\n')
	fw.write(str(n_b) + '\n')
	fw.write(str(dt) + '\n')
	if _image:
		fw.write(str(1) + '\n')
	else:
		fw.write(str(0) + '\n')
	
	if _uniform:
		fw.write(str(1) + '\n')
		fw.write(str(Ub) + '\n')
	elif _shear:
		fw.write(str(2) + '\n')
		fw.write(str(maxV) + '\n')
		fw.write(str(yShear) + '\n')
	elif _BL:
		fw.write(str(3) + '\n')
		fw.write(str(len(yPlus)) + '\n')
		
		for i in range(len(yPlus)):
			fw.write(str(yPlus[i]) + ' ' + str(uVelBL[i]) + ' ' + str(vVelBL[i]) + ' ' + str(wVelBL[i]) + '\n')
		
	else:
		fw.write(str(0) + '\n')
	
	for i in range(0, NP+2):
		
		fw.write(str(Ux[i, j]) + ' ' + str(Uy[i, j]) + ' ' + str(Uz[i, j]) + ' ' + str(Ux_s[i, j]) + ' ' + str(Uy_s[i, j]) + ' ' + str(Uz_s[i, j]) + ' ' + str(SIGMA[i, j]) + '\n')
		
	fw.close()

# This is the main function. It loops over every node of the filament
# and calculates the velocity. vec_local calculates the velocity due to the
# curvature of the filament itself in the binormal direction. 
# In _step == 0, RK45 is used to estimate velocity initially.
# Later steps use a 5th order Adam-Bashforth method. This combination ensures 
# that the results are obtained fast with small numerical errors.
# All of these steps are directly handled by FORTRAN.
# NOTE: FORTRAN is meant only for the M1 KK method. For LIA, use python
# preferably with numba.
def computeNodes(j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, \
beta, SIGMA, delta_ttm, sigma_1, sigma_2, gamma_param, sigma3_1, sigma3_2, ds, n_b, _dtype, \
Vx_0, Vy_0, Vz_0, Vx_m, Vy_m, Vz_m, Vx_m1, Vy_m1, Vz_m1, Vx_m2, Vy_m2, Vz_m2, dt, _LIA, one_o_2ds, one_o_ds2,\
_uniform, _shear, _BL, Ub, yShear, maxV, yPlus, uVelBL, vVelBL, wVelBL, _image):
	
	if _LIA:
		
		print('set _fullFORTRAN = False and run again..')
		raise SystemError
	
	if not _LIA:
		
		# Write data to file
		
		fortran_writer(NP, delta_ttm, sigma_1, sigma_2, gamma_param, j, ds, n_b, Ux, Uy, \
		Uz, Ux_s, Uy_s, Uz_s, SIGMA, _image, _uniform, Ub, _shear, maxV, yShear, _BL, yPlus,\
		uVelBL, vVelBL, wVelBL, dt)
		
		# raise SystemError
		
		os.chdir('src/full_fortran')
		os.system('./uttm_full')

