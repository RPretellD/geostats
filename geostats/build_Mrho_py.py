""" A suite of functions to compute correlation matrices for geostatistical analysis.
"""

import numpy as np
from vincenty import vincenty

__author__ = 'A. Renmin Pretell Ductram'

#===================================================================================================
# rho_E: Sampled location to sampled location
#===================================================================================================
def MrhoE_sam2sam(data_lat, data_lon, L_E, gamma_E):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	
	Returns
	=======
	rho_E correlation matrix.
	"""

	data_lat_arr = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr = np.asarray(data_lon, dtype=np.float64)

	M      = data_lat_arr.size
	MrhoE  = np.empty([M,M], dtype='float64')
	np.fill_diagonal(MrhoE, 1.0)

	inv_LE = 1/L_E

	for i in range(M):
		for j in range(i+1,M):
			d_E        = vincenty((data_lat_arr[i],data_lon_arr[i]),(data_lat_arr[j],data_lon_arr[j]))
			rho_E      = np.exp(-(d_E*inv_LE)**gamma_E)
			MrhoE[i,j] = rho_E
			MrhoE[j,i] = rho_E
	return MrhoE

#===================================================================================================
# rho_EA: Sampled location to sampled location
#===================================================================================================
def MrhoEA_sam2sam(data_lat, data_lon, epi_lat, epi_lon, L_E, gamma_E, L_A):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	epi_lat/epi_lon: Epicenter coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length (deg).
	
	Returns
	=======
	rho_EA correlation matrix.
	"""

	data_lat_arr = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr = np.asarray(data_lon, dtype=np.float64)

	data_lat_arr_r = np.deg2rad(data_lat_arr)
	data_lon_arr_r = np.deg2rad(data_lon_arr)

	epi_lat = np.deg2rad(epi_lat)
	epi_lon = np.deg2rad(epi_lon)

	M      = data_lat_arr.size
	MrhoEA = np.empty([M,M], dtype='float64')
	np.fill_diagonal(MrhoEA, 1.0)

	inv_LE = 1/L_E
	inv_LA = 1/L_A

	for i in range(M):
		lat1 = data_lat_arr_r[i]
		lon1 = data_lon_arr_r[i]
		Az1  = np.arctan2(np.sin(lon1-epi_lon)*np.cos(lat1),np.cos(epi_lat)*np.sin(lat1)-np.sin(epi_lat)*np.cos(lat1)*np.cos(lon1-epi_lon))
		for j in range(i+1,M):
			d_E         = vincenty((data_lat_arr[i],data_lon_arr[i]),(data_lat_arr[j],data_lon_arr[j]))
			rho_E       = np.exp(-(d_E*inv_LE)**gamma_E)
			lat2        = data_lat_arr_r[j]
			lon2        = data_lon_arr_r[j]
			Az2         = np.arctan2(np.sin(lon2-epi_lon)*np.cos(lat2),np.cos(epi_lat)*np.sin(lat2)-np.sin(epi_lat)*np.cos(lat2)*np.cos(lon2-epi_lon))
			d_A         = np.rad2deg(np.arccos(np.cos(Az1-Az2)))
			rho_A       = (1+d_A*inv_LA)*((1-d_A/180)**(180*inv_LA))
			rho_EA      = rho_E*rho_A
			MrhoEA[i,j] = rho_EA
			MrhoEA[j,i] = rho_EA
	return MrhoEA

#===================================================================================================
# rho_EAS: Sampled location to sampled location
#===================================================================================================
def MrhoEAS_sam2sam(data_lat, data_lon, data_Vs30, epi_lat, epi_lon, L_E, gamma_E, L_A, L_S, weight):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	data_Vs30: Sampled location Vs30 (m/s).
	epi_lat/epi_lon: Epicenter coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length (deg).
	L_S: Vs30 dissimilarity distance correlation length (m/s).
	weight: Weigth parameter.
	
	Returns
	=======
	rho_EAS correlation matrix.
	"""

	data_lat_arr  = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr  = np.asarray(data_lon, dtype=np.float64)
	data_Vs30_arr = np.asarray(data_Vs30, dtype=np.float64)

	data_lat_arr_r = np.deg2rad(data_lat_arr)
	data_lon_arr_r = np.deg2rad(data_lon_arr)

	epi_lat = np.deg2rad(epi_lat)
	epi_lon = np.deg2rad(epi_lon)

	M       = data_lat_arr.size
	MrhoEAS = np.empty([M,M], dtype='float64')
	np.fill_diagonal(MrhoEAS, 1.0)

	inv_LE = 1/L_E
	inv_LA = 1/L_A
	inv_LS = 1/L_S

	for i in range(M):
		lat1 = data_lat_arr_r[i]
		lon1 = data_lon_arr_r[i]
		Az1  = np.arctan2(np.sin(lon1-epi_lon)*np.cos(lat1),np.cos(epi_lat)*np.sin(lat1)-np.sin(epi_lat)*np.cos(lat1)*np.cos(lon1-epi_lon))
		for j in range(i+1,M):
			d_E          = vincenty((data_lat_arr[i],data_lon_arr[i]),(data_lat_arr[j],data_lon_arr[j]))
			rho_E        = np.exp(-(d_E*inv_LE)**gamma_E)
			lat2         = data_lat_arr_r[j]
			lon2         = data_lon_arr_r[j]
			Az2          = np.arctan2(np.sin(lon2-epi_lon)*np.cos(lat2),np.cos(epi_lat)*np.sin(lat2)-np.sin(epi_lat)*np.cos(lat2)*np.cos(lon2-epi_lon))
			d_A          = np.rad2deg(np.arccos(np.cos(Az1-Az2)))
			rho_A        = (1+d_A*inv_LA)*((1-d_A/180)**(180*inv_LA))
			d_S          = abs(data_Vs30_arr[i]-data_Vs30_arr[j])
			rho_S        = np.exp(-(d_S*inv_LS))
			rho_EAS      = rho_E*(weight*rho_A+(1-weight)*rho_S)
			MrhoEAS[i,j] = rho_EAS
			MrhoEAS[j,i] = rho_EAS
	return MrhoEAS

#===================================================================================================
# rho_E: Sampled location to unsampled location
#===================================================================================================
def MrhoE_sam2uns(data_lat, data_lon, unkn_lat, unkn_lon, L_E, gamma_E):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	unkn_lat/unkn_lon: Unsampled location coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	
	Returns
	=======
	rhoE correlation matrix.
	"""

	data_lat_arr = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr = np.asarray(data_lon, dtype=np.float64)

	unkn_lat_arr = np.atleast_1d(unkn_lat)
	unkn_lon_arr = np.atleast_1d(unkn_lon)

	M = data_lat_arr.size
	N = unkn_lat_arr.size
	MrhoE = np.empty([M,N], dtype='float64')

	inv_LE = 1/L_E

	for i in range(M):
		for j in range(N):
			d_E        = vincenty((data_lat_arr[i],data_lon_arr[i]),(unkn_lat_arr[j],unkn_lon_arr[j]))
			rhoE       = np.exp(-(d_E*inv_LE)**gamma_E)
			MrhoE[i,j] = rhoE
	return MrhoE

#===================================================================================================
# rho_EA: Sampled location to unsampled location
#===================================================================================================
def MrhoEA_sam2uns(data_lat, data_lon, epi_lat, epi_lon, unkn_lat, unkn_lon, L_E, gamma_E, L_A):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	unkn_lat/unkn_lon: Unsampled location coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length (deg).
	epi_lat/epi_lon: Epicenter coordinates (deg).
	
	Returns
	=======
	rho_EA correlation matrix.
	"""

	data_lat_arr = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr = np.asarray(data_lon, dtype=np.float64)

	unkn_lat_arr = np.atleast_1d(unkn_lat)
	unkn_lon_arr = np.atleast_1d(unkn_lon)

	data_lat_arr_r = np.deg2rad(data_lat_arr)
	data_lon_arr_r = np.deg2rad(data_lon_arr)

	unkn_lat_arr_r = np.deg2rad(unkn_lat_arr)
	unkn_lon_arr_r = np.deg2rad(unkn_lon_arr)

	epi_lat = np.deg2rad(epi_lat)
	epi_lon = np.deg2rad(epi_lon)

	M = data_lat_arr.size
	N = unkn_lat_arr.size
	MrhoEA = np.empty([M,N], dtype='float64')

	inv_LE = 1/L_E
	inv_LA = 1/L_A

	for i in range(N):
		lat1 = unkn_lat_arr_r[i]
		lon1 = unkn_lon_arr_r[i]
		Az1  = np.arctan2(np.sin(lon1-epi_lon)*np.cos(lat1),np.cos(epi_lat)*np.sin(lat1)-np.sin(epi_lat)*np.cos(lat1)*np.cos(lon1-epi_lon))
		for j in range(M):
			d_E         = vincenty((data_lat_arr[j],data_lon_arr[j]),(unkn_lat_arr[i],unkn_lon_arr[i]))
			rho_E       = np.exp(-(d_E*inv_LE)**gamma_E)
			lat2        = data_lat_arr_r[j]
			lon2        = data_lon_arr_r[j]
			Az2         = np.arctan2(np.sin(lon2-epi_lon)*np.cos(lat2),np.cos(epi_lat)*np.sin(lat2)-np.sin(epi_lat)*np.cos(lat2)*np.cos(lon2-epi_lon))
			d_A         = np.rad2deg(np.arccos(np.cos(Az1-Az2)))
			rho_A       = (1+d_A*inv_LA)*((1-d_A/180)**(180*inv_LA))
			rho_EA      = rho_E*rho_A
			MrhoEA[j,i] = rho_EA
	return MrhoEA

#===================================================================================================
# rho_EAS: Sampled location to unsampled location
#===================================================================================================
def MrhoEAS_sam2uns(data_lat, data_lon, data_Vs30, epi_lat, epi_lon, unkn_lat, unkn_lon, unkn_Vs30, L_E, gamma_E, L_A, L_S, weight):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	data_Vs30: Sampled location Vs30 (m/s).
	epi_lat/epi_lon: Epicenter coordinates in deg.
	unkn_lat/unkn_lon: Unsampled location coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length (deg).
	L_S: Vs30 dissimilarity distance correlation length (m/s).
	weight: Weight model coefficient.
	
	Returns
	=======
	rho_EAS correlation matrix.
	"""

	data_lat_arr  = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr  = np.asarray(data_lon, dtype=np.float64)
	data_Vs30_arr = np.asarray(data_Vs30, dtype=np.float64)

	unkn_lat_arr  = np.atleast_1d(unkn_lat)
	unkn_lon_arr  = np.atleast_1d(unkn_lon)
	unkn_Vs30_arr = np.atleast_1d(unkn_Vs30)

	data_lat_arr_r = np.deg2rad(data_lat_arr)
	data_lon_arr_r = np.deg2rad(data_lon_arr)

	unkn_lat_arr_r = np.deg2rad(unkn_lat_arr)
	unkn_lon_arr_r = np.deg2rad(unkn_lon_arr)

	epi_lat = np.deg2rad(epi_lat)
	epi_lon = np.deg2rad(epi_lon)

	M = data_lat_arr.size
	N = unkn_lat_arr.size
	MrhoEAS = np.empty([M,N], dtype='float64')

	inv_LE = 1/L_E
	inv_LA = 1/L_A
	inv_LS = 1/L_S

	for i in range(N):
		lat1 = unkn_lat_arr_r[i]
		lon1 = unkn_lon_arr_r[i]
		Az1  = np.arctan2(np.sin(lon1-epi_lon)*np.cos(lat1),np.cos(epi_lat)*np.sin(lat1)-np.sin(epi_lat)*np.cos(lat1)*np.cos(lon1-epi_lon))
		for j in range(M):
			d_E          = vincenty((data_lat_arr[j],data_lon_arr[j]),(unkn_lat_arr[i],unkn_lon_arr[i]))
			rho_E        = np.exp(-(d_E*inv_LE)**gamma_E)
			lat2         = data_lat_arr_r[j]
			lon2         = data_lon_arr_r[j]
			Az2          = np.arctan2(np.sin(lon2-epi_lon)*np.cos(lat2),np.cos(epi_lat)*np.sin(lat2)-np.sin(epi_lat)*np.cos(lat2)*np.cos(lon2-epi_lon))
			d_A          = np.rad2deg(np.arccos(np.cos(Az1-Az2)))
			rho_A        = (1+d_A*inv_LA)*((1-d_A/180)**(180*inv_LA))
			d_S          = abs(unkn_Vs30_arr[i]-data_Vs30_arr[j])
			rho_S        = np.exp(-(d_S*inv_LS))
			rho_EAS      = rho_E*(weight*rho_A+(1-weight)*rho_S)
			MrhoEAS[j,i] = rho_EAS
	return MrhoEAS

#===================================================================================================
# rho_E: Sampled to unsampled grid
#===================================================================================================
def MrhoE_sam2grid(data_lat, data_lon, grid_lat, grid_lon, L_E, gamma_E):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	grid_lat/grid_lon: Unsampled grid location coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	
	Returns
	=======
	rho_E correlation matrix.
	"""

	data_lat_arr = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr = np.asarray(data_lon, dtype=np.float64)

	grid_lat_arr = np.asarray(grid_lat, dtype=np.float64)
	grid_lon_arr = np.asarray(grid_lon, dtype=np.float64)

	M      = data_lat_arr.size
	n_lats = grid_lat_arr.size
	n_lons = grid_lon_arr.size
	MrhoE  = np.empty([M,n_lats*n_lons], dtype='float64')

	inv_LE = 1/L_E

	for i in range(n_lons):
		for j in range(n_lats):
			for k in range(M):
				d_E                 = vincenty((data_lat_arr[k],data_lon_arr[k]),(grid_lat_arr[j],grid_lon_arr[i]))
				rho_E               = np.exp(-(d_E*inv_LE)**gamma_E)
				MrhoE[k,i*n_lats+j] = rho_E
	return MrhoE

#===================================================================================================
# rho_EA: Sampled to unsampled grid
#===================================================================================================
def MrhoEA_sam2grid(data_lat, data_lon, epi_lat, epi_lon, grid_lat, grid_lon, L_E, gamma_E, L_A):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates (deg).
	grid_lat/grid_lon: Unsampled grid location coordinates (deg).
	L_E: Euclidean distance correlation length (km).
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length (deg).
	epi_lat/epi_lon: Epicenter coordinates (deg).
	
	Returns
	=======
	rho_EA correlation matrix.
	"""

	data_lat_arr = np.asarray(data_lat, dtype=np.float64)
	data_lon_arr = np.asarray(data_lon, dtype=np.float64)

	grid_lat_arr = np.asarray(grid_lat, dtype=np.float64)
	grid_lon_arr = np.asarray(grid_lon, dtype=np.float64)

	epi_lat_r = np.deg2rad(epi_lat)
	epi_lon_r = np.deg2rad(epi_lon)

	grid_lat_r = np.deg2rad(grid_lat_arr)
	grid_lon_r = np.deg2rad(grid_lon_arr)
	data_lat_r = np.deg2rad(data_lat_arr)
	data_lon_r = np.deg2rad(data_lon_arr)

	M      = data_lat_arr.size
	n_lats = grid_lat_arr.size
	n_lons = grid_lon_arr.size
	MrhoEA = np.empty([M,n_lats*n_lons], dtype='float64')

	inv_LE = 1/L_E
	inv_LA = 1/L_A

	for i in range(n_lons):
		for j in range(n_lats):
			lon1 = grid_lon_r[i]
			lat1 = grid_lat_r[j]
			Az1  = np.arctan2(np.sin(lon1-epi_lon_r)*np.cos(lat1),np.cos(epi_lat_r)*np.sin(lat1)-np.sin(epi_lat_r)*np.cos(lat1)*np.cos(lon1-epi_lon_r))
			for k in range(M):
				d_E                  = vincenty((data_lat_arr[k],data_lon_arr[k]),(grid_lat_arr[j],grid_lon_arr[i]))
				rho_E                = np.exp(-(d_E*inv_LE)**gamma_E)
				lat2                 = data_lat_r[k]
				lon2                 = data_lon_r[k]
				Az2                  = np.arctan2(np.sin(lon2-epi_lon_r)*np.cos(lat2),np.cos(epi_lat_r)*np.sin(lat2)-np.sin(epi_lat_r)*np.cos(lat2)*np.cos(lon2-epi_lon_r))
				d_A                  = np.rad2deg(np.arccos(np.cos(Az1-Az2)))
				rho_A                = (1+d_A*inv_LA)*((1-d_A/180)**(180*inv_LA))
				rho_EA               = rho_E*rho_A
				MrhoEA[k,i*n_lats+j] = rho_EA
	return MrhoEA
