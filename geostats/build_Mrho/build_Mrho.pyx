""" A suite of functions to compute correlation matrices for geostatistical analysis.
"""

import numpy as np
cimport numpy as np
cimport cython
from ..geostats_tools.geostats_tools import vincenty_cy
from libc.math cimport pi, sin, cos, tan, atan, atan2, acos, pow, sqrt, exp

__author__ = 'A. Renmin Pretell Ductram'

#===================================================================================================
# rho_E: Sampled location to sampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoE_sam2sam_cy(double[:] data_lat,double[:] data_lon,double L_E,double gamma_E):
	
	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	
	Returns
	=======
	rho_E correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef double [:, :] MrhoE = np.zeros([M,M], dtype='float64')
	cdef int i, j
	
	for i in range(M):
		for j in range(i,M,1):
			d_E   = vincenty_cy(data_lat[i],data_lon[i],data_lat[j],data_lon[j])
			rhoE = exp(-pow(d_E/L_E,gamma_E))
			MrhoE[i][j] = rhoE
			MrhoE[j][i] = MrhoE[i][j]
	return np.asarray(MrhoE)

#===================================================================================================
# rho_EA: Sampled location to sampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEA_sam2sam_cy(double[:] data_lat,double[:] data_lon,double epi_lat,double epi_lon,double L_E,double gamma_E,double L_A):
	
	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	epi_lat/epi_lon: Epicenter coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length in deg.
	
	Returns
	=======
	rho_EA correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef double [:, :] MrhoEA = np.zeros([M,M], dtype='float64')
	cdef int i, j
	cdef double lat1, lat2, lon1, lon2, Az1, Az2, d_E, d_A, rho_E, rho_A
	
	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180
	for i in range(M):
		lat1 = data_lat[i]*pi/180
		lon1 = data_lon[i]*pi/180
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for j in range(i,M,1):
			d_E   = vincenty_cy(data_lat[i],data_lon[i],data_lat[j],data_lon[j])
			rho_E = exp(-pow(d_E/L_E,gamma_E))
			lat2  = data_lat[j]*pi/180
			lon2  = data_lon[j]*pi/180
			Az2   = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A   = acos(cos(Az1-Az2))*180/pi
			rho_A = (1+d_A/L_A)*pow(1-d_A/180,180/L_A)
			MrhoEA[i][j] = rho_E*rho_A
			MrhoEA[j][i] = MrhoEA[i][j]
	return np.asarray(MrhoEA)

#===================================================================================================
# rho_EAS: Sampled location to sampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEAS_sam2sam_cy(double[:] data_lat,double[:] data_lon,double[:] data_Vs30,double epi_lat,double epi_lon,double L_E,double gamma_E,double L_A,double L_S,double weight):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	data_Vs30: Sampled location Vs30.
	epi_lat/epi_lon: Epicenter coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length in deg.
	L_S: Vs30 dissimilarity distance correlation length.
	weight: Weigth parameter.
	
	Returns
	=======
	rho_EAS correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef double [:, :] MrhoEAS = np.zeros([M,M], dtype='float64')
	cdef int i, j
	cdef double lat1, lon1, Az1, d_E, rhoE, lat2, lon2, Az2, d_A, rhoA, d_S, rhoS
	
	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180
	for i in range(M):
		lat1 = data_lat[i]*pi/180
		lon1 = data_lon[i]*pi/180
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for j in range(i,M,1):
			d_E   = vincenty_cy(data_lat[i],data_lon[i],data_lat[j],data_lon[j])
			rhoE = exp(-pow(d_E/L_E,gamma_E))
			lat2 = data_lat[j]*pi/180
			lon2 = data_lon[j]*pi/180
			Az2  = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A  = acos(cos(Az1-Az2))*180/pi
			rhoA = (1+d_A/L_A)*pow(1-d_A/180,180/L_A)
			d_S  = abs(data_Vs30[i]-data_Vs30[j])
			rhoS = exp(-1*d_S/L_S)
			MrhoEAS[i][j] = rhoE*(weight*rhoA+(1-weight)*rhoS)
			MrhoEAS[j][i] = MrhoEAS[i][j]
	return np.asarray(MrhoEAS)

#===================================================================================================
# rho_E: Sampled location to unsampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoE_sam2uns_cy(double[:] data_lat,double[:] data_lon,double[:] unkn_lat,double[:] unkn_lon,double L_E,double gamma_E):
	
	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	unkn_lat/unkn_lon: Unsampled location coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	
	Returns
	=======
	rhoE correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef int N = int(len(unkn_lat))
	cdef double [:, :] MrhoE = np.zeros([M,N], dtype='float64')
	cdef int i, k
	cdef double d_E, rhoE
	
	for i in range(M):
		for k in range(N):
			d_E  = vincenty_cy(data_lat[k],data_lon[k],unkn_lat[i],unkn_lon[i])
			rhoE = exp(-pow(d_E/L_E,gamma_E))
			MrhoE[k][i] = rhoE
	return np.asarray(MrhoE)

#===================================================================================================
# CORRELATION MATRIX: STATION to KRIGING SITE
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEA_sam2uns_cy(double[:] data_lat,double[:] data_lon,double epi_lat,double epi_lon,double[:] unkn_lat,double[:] unkn_lon,double L_E,double gamma_E,double L_A):
	
	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	unkn_lat/unkn_lon: Unsampled location coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length in deg.
	epi_lat/epi_lon: Epicenter coordinates in deg.
	
	Returns
	=======
	rho_EA correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef int N = int(len(unkn_lat))
	cdef double [:, :] MrhoEA = np.zeros([M,N], dtype='float64')
	cdef int i, k
	cdef double lat1, lat2, lon1, lon2, Az1, Az2, d_E, d_A, rhoE, rhoA
	
	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180
	for i in range(N):
		lat1 = unkn_lat[i]*pi/180
		lon1 = unkn_lon[i]*pi/180
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for k in range(M):
			d_E  = vincenty_cy(data_lat[k],data_lon[k],unkn_lat[i],unkn_lon[i])
			rhoE = exp(-pow(d_E/L_E,gamma_E))
			lat2 = data_lat[k]*pi/180
			lon2 = data_lon[k]*pi/180
			Az2  = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A  = acos(cos(Az1-Az2))*180/pi
			rhoA = (1+d_A/L_A)*pow(1-d_A/180,180/L_A)
			MrhoEA[k][i] = rhoE*rhoA
	return np.asarray(MrhoEA)

#===================================================================================================
# rho_EAS: Sampled location to unsampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEAS_sam2uns_cy(double[:] data_lat,double[:] data_lon,double[:] data_Vs30,double epi_lat,double epi_lon,double[:] unkn_lat,double[:] unkn_lon,double[:] unkn_Vs30,double L_E,double gamma_E,double L_A,double L_S,double weight):
	
	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	data_Vs30: Sampled location Vs30.
	epi_lat/epi_lon: Epicenter coordinates in deg.
	unkn_lat/unkn_lon: Unsampled location coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length in deg.
	L_S: Vs30 dissimilarity distance correlation length.
	weight: Weight model coefficient.
	
	Returns
	=======
	rho_EAS correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef int N = int(len(unkn_lat))
	cdef double [:, :] MrhoEAS = np.zeros([M,N], dtype='float64')
	cdef int i, k
	cdef double lat1, lon1, Az1, d_E, rhoE, lat2, lon2, Az2, d_A, rhoA
	
	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180
	for i in range(N):
		lat1 = unkn_lat[i]*pi/180
		lon1 = unkn_lon[i]*pi/180
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for k in range(M):
			d_E  = vincenty_cy(data_lat[k],data_lon[k],unkn_lat[i],unkn_lon[i])
			rhoE = exp(-pow(d_E/L_E,gamma_E))
			lat2 = data_lat[k]*pi/180
			lon2 = data_lon[k]*pi/180
			Az2  = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A  = acos(cos(Az1-Az2))*180/pi
			rhoA = (1+d_A/L_A)*pow(1-d_A/180,180/L_A)
			d_S  = abs(unkn_Vs30[i]-data_Vs30[k])
			rhoS = exp(-1*d_S/L_S)
			MrhoEAS[k][i] = rhoE*(weight*rhoA+(1-weight)*rhoS)
	return np.asarray(MrhoEAS)

#===================================================================================================
# rho_E: Sampled to unsampled grid
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoE_sam2grid_cy(double[:] data_lat,double[:] data_lon,double[:] grid_lat,double[:] grid_lon,double L_E,double gamma_E):
	
	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	grid_lat/grid_lon: Unsampled grid location coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	
	Returns
	=======
	rho_E correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef int n_lats = int(len(grid_lat))
	cdef int n_lons = int(len(grid_lon))
	cdef double [:, :] MrhoE = np.zeros([M,n_lats*n_lons], dtype='float64')
	cdef int i, j, k
	cdef double lat1, lon1, d_E, rhoE, lat2, lon2
	
	for i in range(n_lats):
		for j in range(n_lons):
			for k in range(M):
				d_E   = vincenty_cy(data_lat[k],data_lon[k],grid_lat[i],grid_lon[j])
				rhoE = exp(-pow(d_E/L_E,gamma_E))
				MrhoE[k][j*n_lats+i] = rhoE
	return np.asarray(MrhoE)

#===================================================================================================
# rho_EA: Sampled to unsampled grid
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEA_sam2grid_cy(double[:] data_lat,double[:] data_lon,double epi_lat,double epi_lon,double[:] grid_lat,double[:] grid_lon,double L_E,double gamma_E,double L_A):

	"""
	Parameter
	=========
	data_lat/data_lon: Sampled location coordinates in deg.
	grid_lat/grid_lon: Unsampled grid location coordinates in deg.
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length in deg.
	epi_lat/epi_lon: Epicenter coordinates in deg.
	
	Returns
	=======
	rho_EA correlation matrix.
	"""
	
	cdef int M = int(len(data_lat))
	cdef int n_lats = int(len(grid_lat))
	cdef int n_lons = int(len(grid_lon))
	cdef double [:, :] MrhoEA = np.zeros([M,n_lats*n_lons], dtype='float64')
	cdef int i, j, k
	cdef double lat1, lat2, lon1, lon2, Az1, Az2, d_E, d_A, rhoE, rhoA

	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180
	for i in range(n_lats):
		for j in range(n_lons):
			lat1 = grid_lat[i]*pi/180
			lon1 = grid_lon[j]*pi/180
			Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
			for k in range(M):
				d_E  = vincenty_cy(data_lat[k],data_lon[k],grid_lat[i],grid_lon[j])
				rhoE = exp(-pow(d_E/L_E,gamma_E))
				lat2 = data_lat[k]*pi/180
				lon2 = data_lon[k]*pi/180
				Az2  = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
				d_A  = acos(cos(Az1-Az2))*180/pi
				rhoA = (1+d_A/L_A)*pow(1-d_A/180,180/L_A)
				MrhoEA[k][j*n_lats+i] = rhoE*rhoA
	return np.asarray(MrhoEA)
