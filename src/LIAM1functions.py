import numpy as np
from numba import jit, prange

# C_V is the circumferential component of the 
# core constant which is determined from the 
# similar core solutions of Ting Klein 1991. 
# NOTE: C_W is not calculated since axial flux is
# set to 0.
@jit(nopython = True)
def C_V(delta_bar):
	return (((1 + 0.577 - np.log(2)) / 2) - np.log(delta_bar))

# This function imposes the periodic boundary conditions.
@jit(nopython = True)
def imposeBC(NP, NF, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA):
	
	for j in range(NF):
		
		left_right = np.array([Ux[NP, j] - Ux[1, j], Uy[NP, j] - Uy[1, j], Uz[NP, j] - Uz[1, j]])
		
		Ux[0, j] = Ux[NP-1, j] - left_right[0]
		Uy[0, j] = Uy[NP-1, j] - left_right[1]
		Uz[0, j] = Uz[NP-1, j] - left_right[2]
		
		Ux[NP+1, j] = Ux[2, j] + left_right[0]
		Uy[NP+1, j] = Uy[2, j] + left_right[1]
		Uz[NP+1, j] = Uz[2, j] + left_right[2]
		
		Ux_s[0, j] = Ux_s[NP-1, j]
		Uy_s[0, j] = Uy_s[NP-1, j]
		Uz_s[0, j] = Uz_s[NP-1, j]
		SIGMA[0, j] = SIGMA[NP-1, j]
		
		Ux_s[NP+1, j] = Ux_s[2, j]
		Uy_s[NP+1, j] = Uy_s[2, j]
		Uz_s[NP+1, j] = Uz_s[2, j]
		SIGMA[NP+1, j] = SIGMA[2, j]
		
		Ux_ss[0, j] = Ux_ss[NP-1, j]
		Uy_ss[0, j] = Uy_ss[NP-1, j]
		Uz_ss[0, j] = Uz_ss[NP-1, j]
		
		Ux_ss[NP+1, j] = Ux_ss[2, j]
		Uy_ss[NP+1, j] = Uy_ss[2, j]
		Uz_ss[NP+1, j] = Uz_ss[2, j]
		
		Ux[NP, j] = Ux[1, j] + left_right[0]
		Uy[NP, j] = Uy[1, j] + left_right[1]
		Uz[NP, j] = Uz[1, j] + left_right[2]
		
		Ux_s[NP, j] = Ux_s[1, j]
		Uy_s[NP, j] = Uy_s[1, j]
		Uz_s[NP, j] = Uz_s[1, j]
		SIGMA[NP, j] = SIGMA[1, j]
		
		Ux_ss[NP, j] = Ux_ss[1, j]
		Uy_ss[NP, j] = Uy_ss[1, j]
		Uz_ss[NP, j] = Uz_ss[1, j]
	
	return Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss

# Both computeFiniteDifferencePT and computeFiniteDifference perform
# the same function which is to calculate the finite difference along
# the nodes. In computeFiniteDifferencePT, an additional parameter i needs to 
# be input indicating the location of the node and this is looped over
# automatically in computeFiniteDifference. Two versions of the functions
# are written for convenience. 
@jit(nopython = True)
def computeFiniteDifferencePT(i, j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2, point_c):
	
	point_m = np.array([Ux[i-1, j], Uy[i-1, j], Uz[i-1, j]])
	point_p = np.array([Ux[i+1, j], Uy[i+1, j], Uz[i+1, j]])
	vec_tmp = point_p - point_m
	
	Ux_s[i, j] = vec_tmp[0] * one_o_2ds
	Uy_s[i, j] = vec_tmp[1] * one_o_2ds
	Uz_s[i, j] = vec_tmp[2] * one_o_2ds
	
	SIGMA[i, j] = np.linalg.norm(vec_tmp) * one_o_2ds
	sigma_tmp = SIGMA[i, j]
	
	one_o_sigma3 = 1 / (sigma_tmp**3)
	
	vec_tmp1 = point_p + point_m
	vec_tmp1 = vec_tmp1 - point_c
	vec_tmp1 = vec_tmp1 - point_c
	
	vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
	
	vec_tmp2[0] = vec_tmp2[0] * one_o_2ds * one_o_ds2
	vec_kbx = vec_tmp2[0] * one_o_sigma3
	
	vec_tmp2[1] = vec_tmp2[1] * one_o_2ds * one_o_ds2
	vec_kby = vec_tmp2[1] * one_o_sigma3
	
	vec_tmp2[2] = vec_tmp2[2] * one_o_2ds * one_o_ds2
	vec_kbz = vec_tmp2[2] * one_o_sigma3
	
	Ux_ss[i, j] = vec_kbx
	Uy_ss[i, j] = vec_kby
	Uz_ss[i, j] = vec_kbz
	
	return Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA

@jit(nopython = True)
def computeFiniteDifference(j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2):
	
	for i in range(1, NP+1):
			
		point_m = np.array([Ux[i-1, j], Uy[i-1, j], Uz[i-1, j]])
		point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]])
		point_p = np.array([Ux[i+1, j], Uy[i+1, j], Uz[i+1, j]])
		vec_tmp = point_p - point_m
		
		Ux_s[i, j] = vec_tmp[0] * one_o_2ds
		Uy_s[i, j] = vec_tmp[1] * one_o_2ds
		Uz_s[i, j] = vec_tmp[2] * one_o_2ds
		
		SIGMA[i, j] = np.linalg.norm(vec_tmp) * one_o_2ds
		sigma_tmp = SIGMA[i, j]
		
		one_o_sigma3 = 1 / (sigma_tmp**3)
		
		vec_tmp1 = point_p + point_m
		vec_tmp1 = vec_tmp1 - point_c
		vec_tmp1 = vec_tmp1 - point_c
		
		vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
		
		vec_tmp2[0] = vec_tmp2[0] * one_o_2ds * one_o_ds2
		vec_kbx = vec_tmp2[0] * one_o_sigma3
		
		vec_tmp2[1] = vec_tmp2[1] * one_o_2ds * one_o_ds2
		vec_kby = vec_tmp2[1] * one_o_sigma3
		
		vec_tmp2[2] = vec_tmp2[2] * one_o_2ds * one_o_ds2
		vec_kbz = vec_tmp2[2] * one_o_sigma3
		
		Ux_ss[i, j] = vec_kbx
		Uy_ss[i, j] = vec_kby
		Uz_ss[i, j] = vec_kbz
		
	return Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA

# Implementation of the local induction approximation with Runge Kutta Fehlberg (RK45) method.
@jit(nopython = True)
def LIA_RK45(i, j, NP, one_o_2ds, one_o_ds2, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, _dtype, beta, dt,\
_uniform, _shear, _BL, Ub, yShear, maxV, yPlus, uVelBL, vVelBL, wVelBL):
	
	point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)

	# Compute derivatives with finite difference
	
	Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA = computeFiniteDifferencePT(i, j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2, point_c)
	
	point_c0 = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	vec_kb = np.array([Ux_ss[i, j], Uy_ss[i, j], Uz_ss[i, j]], dtype = _dtype)
	vec_local = np.array([beta * vec_kb[0], beta * vec_kb[1], beta * vec_kb[2]], dtype = _dtype)
	
	if _uniform:
		vec_local[0] = vec_local[0] + Ub
	if _shear:
		if Uy[i, j] > yShear:
			vec_local[0] = vec_local[0] + maxV
		else:
			vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
			
	if _BL:
		vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
		vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
		vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
	
	k1 = dt * vec_local
	
	point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	point_c = point_c + (k1/4)

	# Compute derivatives with finite difference
	
	Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA = computeFiniteDifferencePT(i, j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2, point_c)
	
	point_c0 = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	vec_kb = np.array([Ux_ss[i, j], Uy_ss[i, j], Uz_ss[i, j]], dtype = _dtype)
	vec_local = np.array([beta * vec_kb[0], beta * vec_kb[1], beta * vec_kb[2]], dtype = _dtype)
	
	if _uniform:
		vec_local[0] = vec_local[0] + Ub
	if _shear:
		if Uy[i, j] > yShear:
			vec_local[0] = vec_local[0] + maxV
		else:
			vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
			
	if _BL:
		vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
		vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
		vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
	
	k2 = dt * vec_local
	
	point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	point_c = point_c + (3*k1/32) + (9*k2/32)

	# Compute derivatives with finite difference
	
	Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA = computeFiniteDifferencePT(i, j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2, point_c)
	
	point_c0 = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	vec_kb = np.array([Ux_ss[i, j], Uy_ss[i, j], Uz_ss[i, j]], dtype = _dtype)
	vec_local = np.array([beta * vec_kb[0], beta * vec_kb[1], beta * vec_kb[2]], dtype = _dtype)
	
	if _uniform:
		vec_local[0] = vec_local[0] + Ub
	if _shear:
		if Uy[i, j] > yShear:
			vec_local[0] = vec_local[0] + maxV
		else:
			vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
			
	if _BL:
		vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
		vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
		vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
	
	k3 = dt * vec_local
	
	point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	point_c = point_c + (1932*k1/2197) - (7200*k2/2197) + (7296*k3/2197)

	# Compute derivatives with finite difference
	
	Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA = computeFiniteDifferencePT(i, j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2, point_c)
	
	point_c0 = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	vec_kb = np.array([Ux_ss[i, j], Uy_ss[i, j], Uz_ss[i, j]], dtype = _dtype)
	vec_local = np.array([beta * vec_kb[0], beta * vec_kb[1], beta * vec_kb[2]], dtype = _dtype)
	
	if _uniform:
		vec_local[0] = vec_local[0] + Ub
	if _shear:
		if Uy[i, j] > yShear:
			vec_local[0] = vec_local[0] + maxV
		else:
			vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
			
	if _BL:
		vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
		vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
		vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
	
	k4 = dt * vec_local
	
	point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	point_c = point_c + (439*k1/216) - (8*k2) + (3680*k3/513) - (845*k4/4104)

	# Compute derivatives with finite difference
	
	Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA = computeFiniteDifferencePT(i, j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2, point_c)
	
	point_c0 = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	vec_kb = np.array([Ux_ss[i, j], Uy_ss[i, j], Uz_ss[i, j]], dtype = _dtype)
	vec_local = np.array([beta * vec_kb[0], beta * vec_kb[1], beta * vec_kb[2]], dtype = _dtype)
	
	if _uniform:
		vec_local[0] = vec_local[0] + Ub
	if _shear:
		if Uy[i, j] > yShear:
			vec_local[0] = vec_local[0] + maxV
		else:
			vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
			
	if _BL:
		vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
		vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
		vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
	
	k5 = dt * vec_local
	
	point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	point_c = point_c - (8*k1/27) + (2*k2) - (3544*k3/2565) + (1859*k4/4104) - (11*k5/40)

	# Compute derivatives with finite difference
	
	Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA = computeFiniteDifferencePT(i, j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2, point_c)
	
	point_c0 = np.array([Ux[i, j], Uy[i, j], Uz[i, j]], dtype = _dtype)
	vec_kb = np.array([Ux_ss[i, j], Uy_ss[i, j], Uz_ss[i, j]], dtype = _dtype)
	vec_local = np.array([beta * vec_kb[0], beta * vec_kb[1], beta * vec_kb[2]], dtype = _dtype)
	
	if _uniform:
		vec_local[0] = vec_local[0] + Ub
	if _shear:
		if Uy[i, j] > yShear:
			vec_local[0] = vec_local[0] + maxV
		else:
			vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
			
	if _BL:
		vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
		vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
		vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
	
	k6 = dt * vec_local
	
	return 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55

# This function calculates the thin tube velocity. 
# It includes the periodic boundary condition treatment 
# described in appendix A of Klein Knio 1995.
# The central part of the domain corresponds to u_1 (eq. A4 in KK95)
# where smoothing is applied and left and right images correspond to 
# u_2 (eq. A5 in KK95). The smoothing effect is neglected from the
# images since L >> \delta. 
# This function is mainly meant for the first estimation of velocity with
# RK45.
@jit(nopython = True)
def uttm(point_c, i, j, NP, Ux, Uy, Uz, SIGMA, Ux_s, Uy_s, Uz_s, delta_ttm, \
sigma_1, sigma_2, sigma3_1, sigma3_2, gamma_param, ds, _dtype, n_b):
	
	vec_auto = np.array([0.0, 0.0, 0.0])
	
	ux_tmp = np.zeros((NP+1), dtype = _dtype)
	uy_tmp = np.zeros((NP+1), dtype = _dtype)
	uz_tmp = np.zeros((NP+1), dtype = _dtype)
	
	ux_s = np.zeros((NP+1), dtype = _dtype)
	uy_s = np.zeros((NP+1), dtype = _dtype)
	uz_s = np.zeros((NP+1), dtype = _dtype)
	
	sigma_u = np.zeros((NP+1), dtype = _dtype)
	
	left_right = np.array([Ux[NP, j] - Ux[1, j], Uy[NP, j] - Uy[1, j], Uz[NP, j] - Uz[1, j]])
	flag_overflow = 0
	
	for i_j in range(1, int((NP+1)/2)+1):
		ll = i_j + i - 1
		if ll > NP:
			ll = ll - (NP-1)
			flag_overflow = 1
		if flag_overflow == 0:
			ux_tmp[int(((NP+1)/2)+i_j-1)] = Ux[ll, j]
			uy_tmp[int(((NP+1)/2)+i_j-1)] = Uy[ll, j]
			uz_tmp[int(((NP+1)/2)+i_j-1)] = Uz[ll, j]
		else:
			ux_tmp[int(((NP+1)/2)+i_j-1)] = Ux[ll, j] + left_right[0]
			uy_tmp[int(((NP+1)/2)+i_j-1)] = Uy[ll, j] + left_right[1]
			uz_tmp[int(((NP+1)/2)+i_j-1)] = Uz[ll, j] + left_right[2]
		ux_s[int(((NP+1)/2)+i_j-1)] = Ux_s[ll, j]
		uy_s[int(((NP+1)/2)+i_j-1)] = Uy_s[ll, j]
		uz_s[int(((NP+1)/2)+i_j-1)] = Uz_s[ll, j]
		sigma_u[int(((NP+1)/2)+i_j-1)] = SIGMA[ll, j]
		
	
	flag_overflow = 0
	
	for i_j in range(1, int(((NP+1)/2)-1)+1):
		ll = i - i_j
		if ll < 1:
			ll = (NP-1) + ll
			flag_overflow = 1
		if flag_overflow == 0:
			ux_tmp[int(((NP+1)/2)-i_j)] = Ux[ll, j]
			uy_tmp[int(((NP+1)/2)-i_j)] = Uy[ll, j]
			uz_tmp[int(((NP+1)/2)-i_j)] = Uz[ll, j]
		else:
			ux_tmp[int(((NP+1)/2)-i_j)] = Ux[ll, j] - left_right[0]
			uy_tmp[int(((NP+1)/2)-i_j)] = Uy[ll, j] - left_right[1]
			uz_tmp[int(((NP+1)/2)-i_j)] = Uz[ll, j] - left_right[2]
		ux_s[int(((NP+1)/2)-i_j)] = Ux_s[ll, j]
		uy_s[int(((NP+1)/2)-i_j)] = Uy_s[ll, j]
		uz_s[int(((NP+1)/2)-i_j)] = Uz_s[ll, j]
		sigma_u[int(((NP+1)/2)-i_j)] = SIGMA[ll, j]
		
	delta_ttm_tmp = delta_ttm
	ln_o_ln_coeff = np.log(sigma_1/delta_ttm_tmp)/np.log(sigma_2/sigma_1)
	gamma_tmp = gamma_param[j] / (4 * np.pi)
	
	left_right = np.array([ux_tmp[NP] - ux_tmp[1], uy_tmp[NP] - uy_tmp[1], uz_tmp[NP] - uz_tmp[1]])
	
	# Central part of domain
	
	for i_j in range(1, NP+1):
		
		if not i_j == ((NP+1)/2):
			
			vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
			vec_tmp1 = np.array([ux_tmp[i_j], uy_tmp[i_j], uz_tmp[i_j]])
			
			vec_tmp1 = point_c - vec_tmp1
			one_o_norm3 = np.linalg.norm(vec_tmp1)
			norm3 = one_o_norm3**3
			one_o_norm3 = 1 / one_o_norm3**3
			kernel_tanh_coeff_1 = np.tanh(norm3/sigma3_1)
			kernel_tanh_coeff_2 = np.tanh(norm3/sigma3_2)
			
			vec_tmp1 = one_o_norm3 * vec_tmp1
			
			vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
			vec_auto1 = kernel_tanh_coeff_1 * vec_tmp2
			vec_auto2 = kernel_tanh_coeff_2 * vec_tmp2
			
			vec_tmp2 = vec_auto1 + (vec_auto1 - vec_auto2) * ln_o_ln_coeff
			
			vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
			vec_auto[1] += gamma_tmp * ds * vec_tmp2[1]
			vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
		
	
	for i_b in range(1, n_b):
		
		left_right_b = i_b * left_right
		
		# auto induction of left part
		
		for i_j in range(1, NP):
			
			vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
			vec_tmp1 = np.array([ux_tmp[i_j], uy_tmp[i_j], uz_tmp[i_j]]) - left_right_b
			
			vec_tmp1 = point_c - vec_tmp1
			one_o_norm3 = np.linalg.norm(vec_tmp1)
			one_o_norm3 = 1 / one_o_norm3**3
			
			vec_tmp1 = one_o_norm3 * vec_tmp1
			vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
			
			vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
			vec_auto[1] += gamma_tmp * ds * vec_tmp2[1]
			vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
			
		# auto induction of the right part
		
		for i_j in range(2, NP+1):
			
			vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
			vec_tmp1 = np.array([ux_tmp[i_j], uy_tmp[i_j], uz_tmp[i_j]]) + left_right_b
			
			vec_tmp1 = point_c - vec_tmp1
			one_o_norm3 = np.linalg.norm(vec_tmp1)
			one_o_norm3 = 1 / one_o_norm3**3
			
			vec_tmp1 = one_o_norm3 * vec_tmp1
			vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
			
			vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
			vec_auto[1] += gamma_tmp * ds * vec_tmp2[1]
			vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
	
	if _image:
		# Image part of the vortex
		# Flipping uy_tmp[i_j] gives the image contribution
		# y velocity needs to be subtracted
		
		# Central part of image
		# Even for central part, the effect of smoothing is neglected.
		
		for i_j in range(1, NP+1):
			
			if not i_j == ((NP+1)/2):
				
				vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
				vec_tmp1 = np.array([ux_tmp[i_j], -uy_tmp[i_j], uz_tmp[i_j]])
				
				vec_tmp1 = point_c - vec_tmp1
				one_o_norm3 = np.linalg.norm(vec_tmp1)
				one_o_norm3 = 1 / one_o_norm3**3
				
				vec_tmp1 = one_o_norm3 * vec_tmp1
				vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
				
				vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
				vec_auto[1] -= gamma_tmp * ds * vec_tmp2[1]
				vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
				
		
		for i_b in range(1, n_b):
			
			left_right_b = i_b * left_right
			
			# auto induction of left part
			
			for i_j in range(1, NP):
				
				vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
				vec_tmp1 = np.array([ux_tmp[i_j], -uy_tmp[i_j], uz_tmp[i_j]]) - left_right_b
				
				vec_tmp1 = point_c - vec_tmp1
				one_o_norm3 = np.linalg.norm(vec_tmp1)
				one_o_norm3 = 1 / one_o_norm3**3
				
				vec_tmp1 = one_o_norm3 * vec_tmp1
				vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
				
				vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
				vec_auto[1] -= gamma_tmp * ds * vec_tmp2[1]
				vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
				
			# auto induction of the right part
			
			for i_j in range(2, NP+1):
				
				vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
				vec_tmp1 = np.array([ux_tmp[i_j], -uy_tmp[i_j], uz_tmp[i_j]]) + left_right_b
				
				vec_tmp1 = point_c - vec_tmp1
				one_o_norm3 = np.linalg.norm(vec_tmp1)
				one_o_norm3 = 1 / one_o_norm3**3
				
				vec_tmp1 = one_o_norm3 * vec_tmp1
				vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
				
				vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
				vec_auto[1] -= gamma_tmp * ds * vec_tmp2[1]
				vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
	
	return vec_auto

# This is the main function. It loops over every node of the filament
# and calculates the velocity. vec_local calculates the velocity due to the
# curvature of the filament itself in the binormal direction. 
# In _step == 0, RK45 is used to estimate velocity initially.
# Later steps use a 5th order Adam-Bashforth method. This combination ensures 
# that the results are obtained fast with small numerical errors.
@jit(nopython = True)
def computeNodes(j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, beta, SIGMA, delta_ttm, \
sigma_1, sigma_2, gamma_param, sigma3_1, sigma3_2, ds, n_b, _dtype, Vx, Vy, Vz, _step,\
Vx_0, Vy_0, Vz_0, Vx_m, Vy_m, Vz_m, Vx_m1, Vy_m1, Vz_m1, Vx_m2, Vy_m2, Vz_m2, Ux_0, Uy_0, Uz_0, dt, _LIA, one_o_2ds, one_o_ds2,\
_uniform, _shear, _BL, Ub, yShear, maxV, yPlus, uVelBL, vVelBL, wVelBL, _image):
	
	for i in range(1, NP+1):
			
		point_c0 = np.array([Ux[i, j], Uy[i, j], Uz[i, j]])
		point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]])
		vec_kb = np.array([Ux_ss[i, j], Uy_ss[i, j], Uz_ss[i, j]])
		vec_local = np.array([beta * vec_kb[0], beta * vec_kb[1], beta * vec_kb[2]])
		
		vec_inter = np.array([0.0, 0.0, 0.0])
		
		if not _LIA:
			vec_auto = np.array([0.0, 0.0, 0.0])
			
			ux_tmp = np.zeros((NP+1), dtype = _dtype)
			uy_tmp = np.zeros((NP+1), dtype = _dtype)
			uz_tmp = np.zeros((NP+1), dtype = _dtype)
			
			ux_s = np.zeros((NP+1), dtype = _dtype)
			uy_s = np.zeros((NP+1), dtype = _dtype)
			uz_s = np.zeros((NP+1), dtype = _dtype)
			
			sigma_u = np.zeros((NP+1), dtype = _dtype)
			
			left_right = np.array([Ux[NP, j] - Ux[1, j], Uy[NP, j] - Uy[1, j], Uz[NP, j] - Uz[1, j]])
			flag_overflow = 0
			
			for i_j in range(1, int((NP+1)/2)+1):
				ll = i_j + i - 1
				if ll > NP:
					ll = ll - (NP-1)
					flag_overflow = 1
				if flag_overflow == 0:
					ux_tmp[int(((NP+1)/2)+i_j-1)] = Ux[ll, j]
					uy_tmp[int(((NP+1)/2)+i_j-1)] = Uy[ll, j]
					uz_tmp[int(((NP+1)/2)+i_j-1)] = Uz[ll, j]
				else:
					ux_tmp[int(((NP+1)/2)+i_j-1)] = Ux[ll, j] + left_right[0]
					uy_tmp[int(((NP+1)/2)+i_j-1)] = Uy[ll, j] + left_right[1]
					uz_tmp[int(((NP+1)/2)+i_j-1)] = Uz[ll, j] + left_right[2]
				ux_s[int(((NP+1)/2)+i_j-1)] = Ux_s[ll, j]
				uy_s[int(((NP+1)/2)+i_j-1)] = Uy_s[ll, j]
				uz_s[int(((NP+1)/2)+i_j-1)] = Uz_s[ll, j]
				sigma_u[int(((NP+1)/2)+i_j-1)] = SIGMA[ll, j]
				
			
			flag_overflow = 0
			
			for i_j in range(1, int(((NP+1)/2)-1)+1):
				ll = i - i_j
				if ll < 1:
					ll = (NP-1) + ll
					flag_overflow = 1
				if flag_overflow == 0:
					ux_tmp[int(((NP+1)/2)-i_j)] = Ux[ll, j]
					uy_tmp[int(((NP+1)/2)-i_j)] = Uy[ll, j]
					uz_tmp[int(((NP+1)/2)-i_j)] = Uz[ll, j]
				else:
					ux_tmp[int(((NP+1)/2)-i_j)] = Ux[ll, j] - left_right[0]
					uy_tmp[int(((NP+1)/2)-i_j)] = Uy[ll, j] - left_right[1]
					uz_tmp[int(((NP+1)/2)-i_j)] = Uz[ll, j] - left_right[2]
				ux_s[int(((NP+1)/2)-i_j)] = Ux_s[ll, j]
				uy_s[int(((NP+1)/2)-i_j)] = Uy_s[ll, j]
				uz_s[int(((NP+1)/2)-i_j)] = Uz_s[ll, j]
				sigma_u[int(((NP+1)/2)-i_j)] = SIGMA[ll, j]
				
			
			delta_ttm_tmp = delta_ttm
			ln_o_ln_coeff = np.log(sigma_1/delta_ttm_tmp)/np.log(sigma_2/sigma_1)
			gamma_tmp = gamma_param[j] / (4 * np.pi)
			
			left_right = np.array([ux_tmp[NP] - ux_tmp[1], uy_tmp[NP] - uy_tmp[1], uz_tmp[NP] - uz_tmp[1]])
			
			# Central part of domain
			
			for i_j in range(1, NP+1):
				
				if not i_j == ((NP+1)/2):
					
					vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
					vec_tmp1 = np.array([ux_tmp[i_j], uy_tmp[i_j], uz_tmp[i_j]])
					
					vec_tmp1 = point_c - vec_tmp1
					one_o_norm3 = np.linalg.norm(vec_tmp1)
					norm3 = one_o_norm3**3
					one_o_norm3 = 1 / one_o_norm3**3
					kernel_tanh_coeff_1 = np.tanh(norm3/sigma3_1)
					kernel_tanh_coeff_2 = np.tanh(norm3/sigma3_2)
					
					vec_tmp1 = one_o_norm3 * vec_tmp1
					
					vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
					vec_auto1 = kernel_tanh_coeff_1 * vec_tmp2
					vec_auto2 = kernel_tanh_coeff_2 * vec_tmp2
					
					vec_tmp2 = vec_auto1 + (vec_auto1 - vec_auto2) * ln_o_ln_coeff
					
					vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
					vec_auto[1] += gamma_tmp * ds * vec_tmp2[1]
					vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
				
			
			for i_b in range(1, n_b):
				
				left_right_b = i_b * left_right
				
				# auto induction of left part
				
				for i_j in range(1, NP):
					
					vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
					vec_tmp1 = np.array([ux_tmp[i_j], uy_tmp[i_j], uz_tmp[i_j]]) - left_right_b
					
					vec_tmp1 = point_c - vec_tmp1
					one_o_norm3 = np.linalg.norm(vec_tmp1)
					one_o_norm3 = 1 / one_o_norm3**3
					
					vec_tmp1 = one_o_norm3 * vec_tmp1
					vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
					
					vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
					vec_auto[1] += gamma_tmp * ds * vec_tmp2[1]
					vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
					
				# auto induction of the right part
				
				for i_j in range(2, NP+1):
					
					vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
					vec_tmp1 = np.array([ux_tmp[i_j], uy_tmp[i_j], uz_tmp[i_j]]) + left_right_b
					
					vec_tmp1 = point_c - vec_tmp1
					one_o_norm3 = np.linalg.norm(vec_tmp1)
					one_o_norm3 = 1 / one_o_norm3**3
					
					vec_tmp1 = one_o_norm3 * vec_tmp1
					vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
					
					vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
					vec_auto[1] += gamma_tmp * ds * vec_tmp2[1]
					vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
					
					
			if _image:
				# Image part of the vortex
				# Flipping uy_tmp[i_j] gives the image contribution
				# y velocity needs to be subtracted
				
				# Central part of image
				# Even for central part, the effect of smoothing is neglected.
				
				for i_j in range(1, NP+1):
					
					if not i_j == ((NP+1)/2):
						
						vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
						vec_tmp1 = np.array([ux_tmp[i_j], -uy_tmp[i_j], uz_tmp[i_j]])
						
						vec_tmp1 = point_c - vec_tmp1
						one_o_norm3 = np.linalg.norm(vec_tmp1)
						one_o_norm3 = 1 / one_o_norm3**3
						
						vec_tmp1 = one_o_norm3 * vec_tmp1
						vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
						
						vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
						vec_auto[1] -= gamma_tmp * ds * vec_tmp2[1]
						vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
						
				
				for i_b in range(1, n_b):
					
					left_right_b = i_b * left_right
					
					# auto induction of left part
					
					for i_j in range(1, NP):
						
						vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
						vec_tmp1 = np.array([ux_tmp[i_j], -uy_tmp[i_j], uz_tmp[i_j]]) - left_right_b
						
						vec_tmp1 = point_c - vec_tmp1
						one_o_norm3 = np.linalg.norm(vec_tmp1)
						one_o_norm3 = 1 / one_o_norm3**3
						
						vec_tmp1 = one_o_norm3 * vec_tmp1
						vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
						
						vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
						vec_auto[1] -= gamma_tmp * ds * vec_tmp2[1]
						vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
						
					# auto induction of the right part
					
					for i_j in range(2, NP+1):
						
						vec_tmp = np.array([ux_s[i_j], uy_s[i_j], uz_s[i_j]])
						vec_tmp1 = np.array([ux_tmp[i_j], -uy_tmp[i_j], uz_tmp[i_j]]) + left_right_b
						
						vec_tmp1 = point_c - vec_tmp1
						one_o_norm3 = np.linalg.norm(vec_tmp1)
						one_o_norm3 = 1 / one_o_norm3**3
						
						vec_tmp1 = one_o_norm3 * vec_tmp1
						vec_tmp2 = np.cross(vec_tmp, vec_tmp1)
						
						vec_auto[0] += gamma_tmp * ds * vec_tmp2[0]
						vec_auto[1] -= gamma_tmp * ds * vec_tmp2[1]
						vec_auto[2] += gamma_tmp * ds * vec_tmp2[2]
		
		if _LIA:
			vec_add = vec_inter - vec_local
			# fw = open('vel_LIA.txt', 'a')
			# fw.write(str(vec_local[2]) + '\n')
			# fw.close()
			# print(vec_local[2])
		else:
			vec_add = vec_inter - vec_auto
			# fw = open('vel_M1_'+str(int(sigma_2/sigma_1))+'.txt', 'a')
			# fw.write(str(vec_auto[2]) + '\n')
			# fw.close()
			# print(vec_auto[2])
		
		if _uniform:
			vec_add[0] = vec_add[0] + Ub
		if _shear:
			if Uy[i, j] > yShear:
				vec_add[0] = vec_add[0] + maxV
			else:
				vec_add[0] = vec_add[0] + (Uy[i, j] * maxV / yShear)
				
		if _BL:
			vec_add[0] = vec_add[0] + np.interp(Uy[i, j], yPlus, uVelBL)
			vec_add[1] = vec_add[1] + np.interp(Uy[i, j], yPlus, vVelBL)
			vec_add[2] = vec_add[2] + np.interp(Uy[i, j], yPlus, wVelBL)
		Vx[i, j] = vec_add[0]
		Vy[i, j] = vec_add[1]
		Vz[i, j] = vec_add[2]
		
		if _step == 0:
			
			# Use RK45 for first estimation
			
			if _LIA:
			
				vec_add = LIA_RK45(i, j, NP, one_o_2ds, one_o_ds2, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, _dtype, beta, dt,\
				_uniform, _shear, _BL, Ub, yShear, maxV, yPlus, uVelBL, vVelBL, wVelBL)
			
			else:
				
				vec_local = uttm(point_c, i, j, NP, Ux, Uy, Uz, SIGMA, Ux_s, Uy_s, Uz_s, delta_ttm, \
								sigma_1, sigma_2, sigma3_1, sigma3_2, gamma_param, ds, _dtype, n_b)
				if _uniform:
					vec_local[0] = vec_local[0] + Ub
				if _shear:
					if Uy[i, j] > yShear:
						vec_local[0] = vec_local[0] + maxV
					else:
						vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
						
				if _BL:
					vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
					vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
					vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
				k1 = dt*vec_local
				vec_local = uttm(point_c + (1/4*k1), i, j, NP, Ux, Uy, Uz, SIGMA, Ux_s, Uy_s, Uz_s, delta_ttm, \
								sigma_1, sigma_2, sigma3_1, sigma3_2, gamma_param, ds, _dtype, n_b)
				if _uniform:
					vec_local[0] = vec_local[0] + Ub
				if _shear:
					if Uy[i, j] > yShear:
						vec_local[0] = vec_local[0] + maxV
					else:
						vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
						
				if _BL:
					vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
					vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
					vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
				k2 = dt*vec_local
				vec_local = uttm(point_c + (3/32*k1+9/32*k2), i, j, NP, Ux, Uy, Uz, SIGMA, Ux_s, Uy_s, Uz_s, delta_ttm, \
								sigma_1, sigma_2, sigma3_1, sigma3_2, gamma_param, ds, _dtype, n_b)
				if _uniform:
					vec_local[0] = vec_local[0] + Ub
				if _shear:
					if Uy[i, j] > yShear:
						vec_local[0] = vec_local[0] + maxV
					else:
						vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
						
				if _BL:
					vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
					vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
					vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
				k3 = dt*vec_local
				vec_local = uttm(point_c + (1932/2197*k1-7200/2197*k2 +7296/2197*k3), i, j, NP, Ux, Uy, Uz, SIGMA, Ux_s, Uy_s, Uz_s, delta_ttm, \
								sigma_1, sigma_2, sigma3_1, sigma3_2, gamma_param, ds, _dtype, n_b)
				if _uniform:
					vec_local[0] = vec_local[0] + Ub
				if _shear:
					if Uy[i, j] > yShear:
						vec_local[0] = vec_local[0] + maxV
					else:
						vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
						
				if _BL:
					vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
					vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
					vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
				k4 = dt*vec_local
				vec_local = uttm(point_c + (439/216*k1 -8*k2+3680/513*k3-845/4104*k4), i, j, NP, Ux, Uy, Uz, SIGMA, Ux_s, Uy_s, Uz_s, delta_ttm, \
								sigma_1, sigma_2, sigma3_1, sigma3_2, gamma_param, ds, _dtype, n_b)
				if _uniform:
					vec_local[0] = vec_local[0] + Ub
				if _shear:
					if Uy[i, j] > yShear:
						vec_local[0] = vec_local[0] + maxV
					else:
						vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
						
				if _BL:
					vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
					vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
					vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
				k5 = dt*vec_local
				vec_local = uttm(point_c - 8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5, i, j, NP, Ux, Uy, Uz, SIGMA, Ux_s, Uy_s, Uz_s, delta_ttm, \
								sigma_1, sigma_2, sigma3_1, sigma3_2, gamma_param, ds, _dtype, n_b)
				if _uniform:
					vec_local[0] = vec_local[0] + Ub
				if _shear:
					if Uy[i, j] > yShear:
						vec_local[0] = vec_local[0] + maxV
					else:
						vec_local[0] = vec_local[0] + (Uy[i, j] * maxV / yShear)
						
				if _BL:
					vec_local[0] = vec_local[0] + np.interp(Uy[i, j], yPlus, uVelBL)
					vec_local[1] = vec_local[1] + np.interp(Uy[i, j], yPlus, vVelBL)
					vec_local[2] = vec_local[2] + np.interp(Uy[i, j], yPlus, wVelBL)
				k6 = dt*vec_local
				
				vec_add = 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6
			
			vec_tmp = point_c0 + vec_add
			
			Ux_0[i, j] = vec_tmp[0]
			Uy_0[i, j] = vec_tmp[1]
			Uz_0[i, j] = vec_tmp[2]
			
		elif _step == 1:
			
			c1 = 3/2
			c2 = -1/2
			
			vec_add[0] = dt * (c1 * Vx[i, j] + c2 * Vx_0[i, j])
			vec_add[1] = dt * (c1 * Vy[i, j] + c2 * Vy_0[i, j])
			vec_add[2] = dt * (c1 * Vz[i, j] + c2 * Vz_0[i, j])
			
			vec_tmp = point_c0 + vec_add
			
			Ux_0[i, j] = vec_tmp[0]
			Uy_0[i, j] = vec_tmp[1]
			Uz_0[i, j] = vec_tmp[2]
		
		elif _step == 2:
			
			c1 = 23/12
			c2 = -16/12
			c3 = 5/12
			
			vec_add[0] = dt * (c1 * Vx[i, j] + c2 * Vx_0[i, j] + c3 * Vx_m[i, j])
			vec_add[1] = dt * (c1 * Vy[i, j] + c2 * Vy_0[i, j] + c3 * Vy_m[i, j])
			vec_add[2] = dt * (c1 * Vz[i, j] + c2 * Vz_0[i, j] + c3 * Vz_m[i, j])
			
			vec_tmp = point_c0 + vec_add
			
			Ux_0[i, j] = vec_tmp[0]
			Uy_0[i, j] = vec_tmp[1]
			Uz_0[i, j] = vec_tmp[2]
			
		elif _step == 3:
			
			c1 = 55/24
			c2 = -59/24
			c3 = 37/24
			c4 = -9/24
			
			vec_add[0] = dt * (c1 * Vx[i, j] + c2 * Vx_0[i, j] + c3 * Vx_m[i, j] + c4 * Vx_m1[i, j])
			vec_add[1] = dt * (c1 * Vy[i, j] + c2 * Vy_0[i, j] + c3 * Vy_m[i, j] + c4 * Vy_m1[i, j])
			vec_add[2] = dt * (c1 * Vz[i, j] + c2 * Vz_0[i, j] + c3 * Vz_m[i, j] + c4 * Vz_m1[i, j])
			
			vec_tmp = point_c0 + vec_add
			
			Ux_0[i, j] = vec_tmp[0]
			Uy_0[i, j] = vec_tmp[1]
			Uz_0[i, j] = vec_tmp[2]
		
		else:
			
			c1 = 1901/720
			c2 = -2774/720
			c3 = 2616/720
			c4 = -1274/720
			c5 = 251/720
			
			vec_add[0] = dt * (c1 * Vx[i, j] + c2 * Vx_0[i, j] + c3 * Vx_m[i, j] + c4 * Vx_m1[i, j] + c5 * Vx_m2[i, j])
			vec_add[1] = dt * (c1 * Vy[i, j] + c2 * Vy_0[i, j] + c3 * Vy_m[i, j] + c4 * Vy_m1[i, j] + c5 * Vy_m2[i, j])
			vec_add[2] = dt * (c1 * Vz[i, j] + c2 * Vz_0[i, j] + c3 * Vz_m[i, j] + c4 * Vz_m1[i, j] + c5 * Vz_m2[i, j])
			
			vec_tmp = point_c0 + vec_add
			
			Ux_0[i, j] = vec_tmp[0]
			Uy_0[i, j] = vec_tmp[1]
			Uz_0[i, j] = vec_tmp[2]
		
		Vx_0[i, j] = Vx[i, j]
		Vy_0[i, j] = Vy[i, j]
		Vz_0[i, j] = Vz[i, j]
		
		Vx_m[i, j] = Vx_0[i, j]
		Vy_m[i, j] = Vy_0[i, j]
		Vz_m[i, j] = Vz_0[i, j]
		
		Vx_m1[i, j] = Vx_m[i, j]
		Vy_m1[i, j] = Vy_m[i, j]
		Vz_m1[i, j] = Vz_m[i, j]
		
		Vx_m2[i, j] = Vx_m1[i, j]
		Vy_m2[i, j] = Vy_m1[i, j]
		Vz_m2[i, j] = Vz_m1[i, j]
		
	return Ux_0, Uy_0, Uz_0, Vx_0, Vy_0, Vz_0, Vx_m, Vy_m, Vz_m, Vx_m1, Vy_m1, Vz_m1, Vx_m2, Vy_m2, Vz_m2, Vx, Vy, Vz

# This function is useful when velocity outputs are necessary in between time steps.  
# Numba does not allow writing data to file therefore the decorator is removed. 
# NOTE: This function will execute much slower than computeNodes. However, it is useful
# for debugging purposes.
