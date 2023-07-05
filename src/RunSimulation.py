import numpy as np
import datetime, time

def initializeAndRun(NP, NF, _dtype, ts, u_tau, nu, delta_nu, nsteps, epsilon, nu_bar_param, \
A, alpha, beta, xinit, yinit, L, _uniform, Ub, _shear, maxV, yShear, _BL, \
uVelBL, vVelBL, wVelBL, yPlus, gamma_param_1, m_0_param_1, delta_0_bar_param_1, \
K_M1KK, Phi_M1KK, c_ttm, n_b, _restart, filename, _LIA, overlapRequirement, writeEveryNSteps, \
_partialFORTRAN, _fullFORTRAN, _image, logspaced):
	
	if _partialFORTRAN:
		
		from LIAM1functionsPartialFortran import C_V, imposeBC, computeFiniteDifference, computeNodes
	
	elif _fullFORTRAN:
		
		from LIAM1functionsFullFortran import C_V, computeNodes
	
	else:
	
		try:
			from numba import jit, prange
			from LIAM1functions import C_V, imposeBC, computeFiniteDifference, computeNodes
			# raise SystemError (Use this to force python only mode)
		except:
			print('Installation of numba is recommended.')
			print('NOTE: The program will be very slow without numba.')
			from LIAM1functionsNoNumba import C_V, imposeBC, computeFiniteDifference, computeNodes
	
	# Initialize arrays to store data with the selected _dtype
	
	Ux = np.zeros((NP+2, NF), dtype = _dtype)
	Uy = np.zeros((NP+2, NF), dtype = _dtype)
	Uz = np.zeros((NP+2, NF), dtype = _dtype)

	Ux_s = np.zeros((NP+2, NF), dtype = _dtype)
	Uy_s = np.zeros((NP+2, NF), dtype = _dtype)
	Uz_s = np.zeros((NP+2, NF), dtype = _dtype)

	Ux_ss = np.zeros((NP+2, NF), dtype = _dtype)
	Uy_ss = np.zeros((NP+2, NF), dtype = _dtype)
	Uz_ss = np.zeros((NP+2, NF), dtype = _dtype)

	SIGMA = np.zeros((NP+2, NF), dtype = _dtype)

	S_0 = np.zeros((NF), dtype = _dtype)
	delta_0_bar_param = np.zeros((NF), dtype = _dtype)
	m_0_param = np.zeros((NF), dtype = _dtype)
	gamma_param = np.zeros((NF), dtype = _dtype)
	m_flux = np.zeros((NF), dtype = _dtype)
	int_S_o_S0 = np.zeros((NF), dtype = _dtype)
	S = np.zeros((NF), dtype = _dtype)
	c_v = np.zeros((NF), dtype = _dtype)
	c_w = np.zeros((NF), dtype = _dtype)
	delta_bar = np.zeros((NF), dtype = _dtype)

	print('Memory allocation complete..')
	
	# Overwrite existing restart file, if any.
	
	if not _restart:
		fw = open(filename, 'w+')
		fw.close()
		
	# Set all parameters
	
	gamma_param[0] = gamma_param_1
	m_0_param[0] = m_0_param_1 
	delta_0_bar_param[0] = delta_0_bar_param_1 

	# Create filament or read from existing file

	if not _restart:
		if not logspaced:
			vals = np.linspace(-L, L, NP+2)
		else:
			vals = np.logspace(-1, 1, int((NP/2)+1))/10
			vals = vals*L
			tmp = np.linspace(-L, L, NP+2)
			tmp[:int((NP/2)+1)] = -vals[::-1]
			tmp[int((NP/2)+1)+1:] = vals
			vals = tmp
		for i in range(NP+2):
			Ux[i, 0] = xinit + A*np.cos(alpha)*np.exp(-beta*(vals[i])**2)
			Uy[i, 0] = yinit + A*np.sin(alpha)*np.exp(-beta*(vals[i])**2)
			Uz[i, 0] = (vals[i])
			
			# Write initial data to file
			fw = open(filename, 'a')
			fw.write(str(Ux[i, 0]) + ' ' + str(Uy[i, 0]) + ' ' + str(Uz[i, 0]) + '\n')
			fw.close()
			
	else:
		tmp = np.loadtxt('restart_' + filename.split('.')[0] + '.txt')
		timeRestart = np.loadtxt('time_' + filename.split('.')[0] + '.txt')
		for i in range(NP+2):
			Ux[i, 0] = tmp[i, 0]
			Uy[i, 0] = tmp[i, 1]
			Uz[i, 0] = tmp[i, 2]

	print('Filament has been initialized..')

	# Compute derivative and find length

	dt = ts;
	ds = (2 * np.pi)/(NP-1);
	one_o_2ds = 1 / (ds * 2);
	
	if _fullFORTRAN:
		fw = open('fio.dat', 'w+')
		fw.close()
	
	# Store overlap data
	if _restart:
		fwo = open(filename.split('.')[0] + '_overlap.txt', 'a')
	else:
		fwo = open(filename.split('.')[0] + '_overlap.txt', 'w+')

	for j in range(NF):
		
		maxOverlapVal = 0.0
		
		for i in range(1, NP+1):
			point_m = np.array([Ux[i-1, j], Uy[i-1, j], Uz[i-1, j]])
			point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]])
			point_p = np.array([Ux[i+1, j], Uy[i+1, j], Uz[i+1, j]])
			
			vec_tmp = point_p - point_m
			
			Ux_s[i, j] = vec_tmp[0] * one_o_2ds
			Uy_s[i, j] = vec_tmp[1] * one_o_2ds
			Uz_s[i, j] = vec_tmp[2] * one_o_2ds
			
			SIGMA[i, j] = np.linalg.norm(vec_tmp) * one_o_2ds
			
			if np.linalg.norm(point_c-point_m) > maxOverlapVal:
				maxOverlapVal = np.linalg.norm(point_c-point_m)
		
		if maxOverlapVal < (overlapRequirement):
			print('Maximum overlap value:', maxOverlapVal)
			print('Overlap requirement:', overlapRequirement)
			fwo.write(str(maxOverlapVal) + ' ' + str(epsilon) + ' 1 ' + str(overlapRequirement) + ' 1 \n')
			print('Overlap satisfied')
		else:
			print(maxOverlapVal, epsilon/3)
			if maxOverlapVal < (epsilon):
				fwo.write(str(maxOverlapVal) + ' ' + str(epsilon) + ' 1 ' + str(overlapRequirement) + ' 0 \n')
			else:
				fwo.write(str(maxOverlapVal) + ' ' + str(epsilon) + ' 0 ' + str(overlapRequirement) + ' 0 \n')
			print('Overlap not satisfied')
		
		S_tmp = 0
		
		for i in range(1, NP):
			S_tmp += SIGMA[i, j]
		
		S_0[j] = S_tmp * ds
		
		c_v[j] = C_V(delta_0_bar_param[j])
		delta_bar[j] = delta_0_bar_param[j]
		S0_o_S = 1.0
		
		m_flux[j] = m_0_param[j];
		int_S_o_S0[j] = 0;
		S[j] = S_0[j];

	fwo.close()

	# Initialize for step

	one_o_ds2 = 1 / (ds * ds)

	# Checked. Similar to EZ.

	print('Initialization complete..')
	
	if _fullFORTRAN:
		
		fw = open('fio.dat', 'a')
		fw.write(str(writeEveryNSteps) + '\n')
		fw.write(str(S_0[0]) + '\n')
		fw.write(str(delta_0_bar_param[0]) + '\n')
		fw.write(str(nu_bar_param) + '\n')
		fw.write(str(gamma_param[0]) + '\n')
		fw.write(str(epsilon) + '\n')
		fw.write(str(c_ttm) + '\n')
		fw.write(str(nsteps) + '\n')
		fw.close()
	
	# Main loop

	Vx_0 = np.zeros((NP+2, NF), dtype = _dtype)
	Vy_0 = np.zeros((NP+2, NF), dtype = _dtype)
	Vz_0 = np.zeros((NP+2, NF), dtype = _dtype)

	Vx_m = np.zeros((NP+2, NF), dtype = _dtype)
	Vy_m = np.zeros((NP+2, NF), dtype = _dtype)
	Vz_m = np.zeros((NP+2, NF), dtype = _dtype)

	Vx_m1 = np.zeros((NP+2, NF), dtype = _dtype)
	Vy_m1 = np.zeros((NP+2, NF), dtype = _dtype)
	Vz_m1 = np.zeros((NP+2, NF), dtype = _dtype)

	Vx_m2 = np.zeros((NP+2, NF), dtype = _dtype)
	Vy_m2 = np.zeros((NP+2, NF), dtype = _dtype)
	Vz_m2 = np.zeros((NP+2, NF), dtype = _dtype)
	
	
	if _fullFORTRAN:
		
		beta = gamma_param[j]/(4*np.pi)
		
		# M1 method of Knio Klein
		
		delta_ttm = epsilon
		
		delta_ttm = epsilon
		delta_ttm = delta_ttm * np.exp(c_ttm - c_v[j] + 1 - c_w[j])
		
		sigma_max = np.max(SIGMA[:, j])
		sigma_max = ds * sigma_max
		
		sigma_1 = K_M1KK * sigma_max
		sigma_2 = Phi_M1KK * sigma_1
		
		sigma3_1 = sigma_1**3
		sigma3_2 = sigma_2**3
		
		computeNodes(j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, \
		beta, SIGMA, delta_ttm, sigma_1, sigma_2, gamma_param, sigma3_1, sigma3_2, ds, n_b, _dtype, \
		Vx_0, Vy_0, Vz_0, Vx_m, Vy_m, Vz_m, Vx_m1, Vy_m1, Vz_m1, Vx_m2, Vy_m2, Vz_m2, dt, _LIA, one_o_2ds, one_o_ds2,\
		_uniform, _shear, _BL, Ub, yShear, maxV, yPlus, uVelBL, vVelBL, wVelBL, _image)
	
	else:
		
		for _step in range(nsteps):
			
			start_time = time.time()
			
			# Impose boundary conditions
			
			Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss = imposeBC(NP, NF, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA)
			
			# Compute beta
			
			Vx = np.zeros((NP+2, NF), dtype = _dtype)
			Vy = np.zeros((NP+2, NF), dtype = _dtype)
			Vz = np.zeros((NP+2, NF), dtype = _dtype)
			
			Ux_0 = np.zeros((NP+2, NF), dtype = _dtype)
			Uy_0 = np.zeros((NP+2, NF), dtype = _dtype)
			Uz_0 = np.zeros((NP+2, NF), dtype = _dtype)
			
			for j in range(NF):
				beta = gamma_param[j]/(4*np.pi)
				
				if _LIA:
					beta =  beta * (-np.log(epsilon)+np.log(2)-1+c_v[j])
				
				# M1 method of Knio Klein
				
				delta_ttm = epsilon
				
				delta_ttm = epsilon
				delta_ttm = delta_ttm * np.exp(c_ttm - c_v[j] + 1 - c_w[j])
				
				sigma_max = np.max(SIGMA[:, j])
				sigma_max = ds * sigma_max
				
				sigma_1 = K_M1KK * sigma_max
				sigma_2 = Phi_M1KK * sigma_1
				
				sigma3_1 = sigma_1**3
				sigma3_2 = sigma_2**3
				
				# Compute derivatives with finite difference
				
				Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA = computeFiniteDifference(j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA, one_o_2ds, one_o_ds2)
				
				# Impose BC again
				
				Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss = imposeBC(NP, NF, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, SIGMA)
				
				# Loop over elements
				
				Ux_0, Uy_0, Uz_0, Vx_0, Vy_0, Vz_0, Vx_m, Vy_m, Vz_m, Vx_m1, Vy_m1, Vz_m1, Vx_m2, Vy_m2, Vz_m2, Vx, Vy, Vz = computeNodes(j, NP, Ux, Uy, Uz, Ux_s, Uy_s, Uz_s, Ux_ss, Uy_ss, Uz_ss, \
				beta, SIGMA, delta_ttm, sigma_1, sigma_2, gamma_param, sigma3_1, sigma3_2, ds, n_b, _dtype, Vx, Vy, Vz, _step,\
				Vx_0, Vy_0, Vz_0, Vx_m, Vy_m, Vz_m, Vx_m1, Vy_m1, Vz_m1, Vx_m2, Vy_m2, Vz_m2, Ux_0, Uy_0, Uz_0, dt, _LIA, one_o_2ds, one_o_ds2,\
				_uniform, _shear, _BL, Ub, yShear, maxV, yPlus, uVelBL, vVelBL, wVelBL, _image)
				
				for i in range(1, NP+1):
					
					Ux[i, j] = Ux_0[i, j]
					Uy[i, j] = Uy_0[i, j]
					Uz[i, j] = Uz_0[i, j]
				
				if (_step+1)%writeEveryNSteps == 0:
				
					for i in range(0, NP):
						
						fw = open(filename, 'a')
						fw.write(str(Ux[i, j]) + ' ' + str(Uy[i, j]) + ' ' + str(Uz[i, j]) + '\n')
						fw.close()
				
				# Find new length
				
				S_tmp = 0
				
				for i in range(1, NP):
					S_tmp += SIGMA[i, j]
				
				S[j] = S_tmp * ds
				
				# Find the new stretching terms
				
				S0_o_S = S_0[j] / S[j]
				
				int_S_o_S0[j] = int_S_o_S0[j] + dt * (1 / S0_o_S)
				
				# Find new core thickness
				
				delta_bar[j] = delta_0_bar_param[j]**2
				delta_bar[j] += 4 * nu_bar_param * int_S_o_S0[j]
				delta_bar[j] = np.sqrt(delta_bar[j] * S0_o_S)
				
				m_flux[j] = m_0_param[j] * S0_o_S**2
				
				# Find new C_v
				
				c_v[j] = C_V(delta_bar[j])
				
				# Checked. Similar to EZ.
			if not _restart:
				print('time step:', _step+1, '/', nsteps)
				print('time:', (_step+1) * dt)
				print('Time to finish:', datetime.timedelta(seconds=(time.time() - start_time) * (nsteps - _step - 1)))
			else:
				print('time step:', timeRestart+_step+1, '/', timeRestart+nsteps)
				print('time:', (timeRestart+_step+1) * dt)
				print('Time to finish:', datetime.timedelta(seconds=(time.time() - start_time) * (nsteps - _step - 1)))
				
			# For _step > 0, check overlap condition at each node. 
			
			if _step > 0:
				
				fwo = open(filename.split('.')[0] + '_overlap.txt', 'a')
				
				for i in range(1, NP+1):
					point_m = np.array([Ux[i-1, j], Uy[i-1, j], Uz[i-1, j]])
					point_c = np.array([Ux[i, j], Uy[i, j], Uz[i, j]])
					point_p = np.array([Ux[i+1, j], Uy[i+1, j], Uz[i+1, j]])
					
					vec_tmp = point_p - point_m
					
					Ux_s[i, j] = vec_tmp[0] * one_o_2ds
					Uy_s[i, j] = vec_tmp[1] * one_o_2ds
					Uz_s[i, j] = vec_tmp[2] * one_o_2ds
					
					SIGMA[i, j] = np.linalg.norm(vec_tmp) * one_o_2ds
					
					if np.linalg.norm(point_c-point_m) > maxOverlapVal:
						maxOverlapVal = np.linalg.norm(point_c-point_m)
				
				if maxOverlapVal < (epsilon/3):
					fwo.write(str(maxOverlapVal) + ' ' + str(epsilon) + ' 1 ' + str(epsilon/3) + ' 1 \n')
				else:
					if maxOverlapVal < (epsilon):
						fwo.write(str(maxOverlapVal) + ' ' + str(epsilon) + ' 1 ' + str(epsilon/3) + ' 0 \n')
					else:
						fwo.write(str(maxOverlapVal) + ' ' + str(epsilon) + ' 0 ' + str(epsilon/3) + ' 0 \n')

	# Write restart file
	
	if _fullFORTRAN:
		
		# Get data from fout.dat
		allData = np.loadtxt('fout.dat')
		Ux = allData[:, 0]
		Uy = allData[:, 1]
		Uz = allData[:, 2]
		
		# Get last time step data
		Ux = Ux[-NP-2:]
		Uy = Uy[-NP-2:]
		Uz = Uz[-NP-2:]
		
		import os
		os.chdir('../../')
		

	fw = open('restart_' + filename.split('.')[0] + '.txt', 'w+')
	for i in range(0, NP+2):
		if not _fullFORTRAN:
			fw.write(str(Ux[i, j]) + ' ' + str(Uy[i, j]) + ' ' + str(Uz[i, j]) + '\n')
		else:
			fw.write(str(Ux[i]) + ' ' + str(Uy[i]) + ' ' + str(Uz[i]) + '\n')
	fw.close()

	if not _restart:
		fw = open('time_' + filename.split('.')[0] + '.txt', 'w+')
		if not _fullFORTRAN:
			fw.write(str(_step+1) + '\n')
		else:
			fw.write(str(nsteps) + '\n')
		fw.close()
	else:
		fw = open('time_' + filename.split('.')[0] + '.txt', 'w+')
		if not _fullFORTRAN:
			fw.write(str(timeRestart+_step+1) + '\n')
		else:
			fw.write(str(timeRestart+nsteps) + '\n')
		fw.close()
	
	return Ux, Uy, Uz

