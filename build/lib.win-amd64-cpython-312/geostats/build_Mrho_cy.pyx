""" A suite of functions to compute correlation matrices for geostatistical analysis.
"""

import numpy as np
cimport numpy as np
cimport cython
from .geostats_tools_cy import vincenty_cy
from libc.math cimport sin, cos, tan, atan, atan2, acos, pow, sqrt, exp, fabs

__author__ = 'A. Renmin Pretell Ductram'

#===================================================================================================
# rho_E: Sampled location to sampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoE_sam2sam(data_lat, data_lon, double L_E, double gamma_E):

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

	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr = np.ascontiguousarray(data_lon, dtype=np.float64)
	cdef Py_ssize_t M = data_lat_arr.shape[0]
	cdef np.ndarray[np.double_t, ndim=2] MrhoE = np.empty((M, M), dtype=np.float64)
	cdef Py_ssize_t i, j
	cdef double d_E, rho_E

	cdef double inv_LE = 1.0/L_E

	for i in range(M):
		MrhoE[i,i] = 1.0

	for i in range(M):
		for j in range(i+1,M):
			d_E        = vincenty_cy(data_lat_arr[i],data_lon_arr[i],data_lat_arr[j],data_lon_arr[j])
			rho_E      = exp(-pow(d_E*inv_LE,gamma_E))
			MrhoE[i,j] = rho_E
			MrhoE[j,i] = rho_E
	return MrhoE

#===================================================================================================
# rho_EA: Sampled location to sampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEA_sam2sam(data_lat, data_lon, double epi_lat, double epi_lon, double L_E, double gamma_E, double L_A):

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

	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr = np.ascontiguousarray(data_lon, dtype=np.float64)

	cdef Py_ssize_t M = data_lat_arr.shape[0]
	cdef np.ndarray[np.double_t, ndim=2] MrhoEA = np.empty((M, M), dtype=np.float64)
	cdef Py_ssize_t i, j
	cdef double lat1, lat2, lon1, lon2, Az1, Az2, d_E, d_A, rho_E, rho_A, rho_EA

	cdef double pi = 3.141592653589793
	cdef double deg_to_rad = pi/180.0
	cdef double rad_to_deg = 180.0/pi

	epi_lat *= deg_to_rad
	epi_lon *= deg_to_rad

	cdef double inv_LE = 1.0/L_E
	cdef double inv_LA = 1.0/L_A

	for i in range(M):
		MrhoEA[i,i] = 1.0

	for i in range(M):
		lat1 = data_lat_arr[i]*deg_to_rad
		lon1 = data_lon_arr[i]*deg_to_rad
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for j in range(i+1,M):
			d_E         = vincenty_cy(data_lat_arr[i],data_lon_arr[i],data_lat_arr[j],data_lon_arr[j])
			rho_E       = exp(-pow(d_E*inv_LE,gamma_E))
			lat2        = data_lat_arr[j]*deg_to_rad
			lon2        = data_lon_arr[j]*deg_to_rad
			Az2         = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A         = acos(cos(Az1-Az2))*rad_to_deg
			rho_A       = (1+d_A*inv_LA)*pow(1-d_A/180,180*inv_LA)
			rho_EA      = rho_E*rho_A
			MrhoEA[i,j] = rho_EA
			MrhoEA[j,i] = rho_EA
	return MrhoEA

#===================================================================================================
# rho_EAS: Sampled location to sampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEAS_sam2sam(data_lat, data_lon, data_Vs30, double epi_lat, double epi_lon, double L_E, double gamma_E, double L_A, double L_S, double weight):

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

	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr  = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr  = np.ascontiguousarray(data_lon, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_Vs30_arr = np.ascontiguousarray(data_Vs30, dtype=np.float64)

	cdef Py_ssize_t M = data_lat_arr.shape[0]
	cdef np.ndarray[np.double_t, ndim=2] MrhoEAS = np.empty((M, M), dtype=np.float64)
	cdef Py_ssize_t i, j
	cdef double lat1, lon1, Az1, d_E, rho_E, lat2, lon2, Az2, d_A, rho_A, d_S, rho_S, rho_EAS

	cdef double pi = 3.141592653589793
	cdef double deg_to_rad = pi/180.0
	cdef double rad_to_deg = 180.0/pi

	epi_lat *= deg_to_rad
	epi_lon *= deg_to_rad

	cdef double inv_LE = 1.0/L_E
	cdef double inv_LA = 1.0/L_A
	cdef double inv_LS = 1.0/L_S
	cdef double one_minus_w = 1-weight

	for i in range(M):
		MrhoEAS[i,i] = 1.0

	for i in range(M):
		lat1 = data_lat_arr[i]*deg_to_rad
		lon1 = data_lon_arr[i]*deg_to_rad
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for j in range(i+1,M):
			d_E          = vincenty_cy(data_lat_arr[i],data_lon_arr[i],data_lat_arr[j],data_lon_arr[j])
			rho_E        = exp(-pow(d_E*inv_LE,gamma_E))
			lat2         = data_lat_arr[j]*deg_to_rad
			lon2         = data_lon_arr[j]*deg_to_rad
			Az2          = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A          = acos(cos(Az1-Az2))*rad_to_deg
			rho_A        = (1+d_A*inv_LA)*pow(1-d_A/180,180*inv_LA)
			d_S          = fabs(data_Vs30_arr[i]-data_Vs30_arr[j])
			rho_S        = exp(-(d_S*inv_LS))
			rho_EAS      = rho_E*(weight*rho_A+one_minus_w*rho_S)
			MrhoEAS[i,j] = rho_EAS
			MrhoEAS[j,i] = rho_EAS
	return MrhoEAS

#===================================================================================================
# rho_E: Sampled location to unsampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoE_sam2uns(data_lat, data_lon, unkn_lat, unkn_lon,double L_E,double gamma_E):

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

	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr = np.ascontiguousarray(data_lon, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] unkn_lat_arr = np.ascontiguousarray(unkn_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] unkn_lon_arr = np.ascontiguousarray(unkn_lon, dtype=np.float64)

	cdef Py_ssize_t M = data_lat_arr.shape[0]
	cdef Py_ssize_t N = unkn_lat_arr.shape[0]

	cdef np.ndarray[np.double_t, ndim=2] MrhoE = np.empty((M, N), dtype=np.float64)

	cdef Py_ssize_t i, j
	cdef double d_E, rho_E

	cdef double inv_LE = 1.0/L_E

	for i in range(M):
		for j in range(N):
			d_E  = vincenty_cy(data_lat_arr[i],data_lon_arr[i],unkn_lat_arr[j],unkn_lon_arr[j])
			rho_E = exp(-pow(d_E*inv_LE,gamma_E))
			MrhoE[i,j] = rho_E
	return MrhoE

#===================================================================================================
# rho_EA: Sampled location to unsampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEA_sam2uns(data_lat, data_lon, double epi_lat,double epi_lon, unkn_lat, unkn_lon, double L_E, double gamma_E, double L_A):

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

	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr = np.ascontiguousarray(data_lon, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] unkn_lat_arr = np.ascontiguousarray(unkn_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] unkn_lon_arr = np.ascontiguousarray(unkn_lon, dtype=np.float64)

	cdef Py_ssize_t M = data_lat_arr.shape[0]
	cdef Py_ssize_t N = unkn_lat_arr.shape[0]
	
	cdef np.ndarray[np.double_t, ndim=2] MrhoEA = np.empty((M, N), dtype=np.float64)

	cdef Py_ssize_t i, j
	cdef double lat1, lat2, lon1, lon2, Az1, Az2, d_E, d_A, rho_E, rho_A, rho_EA

	cdef double pi = 3.141592653589793
	cdef double deg_to_rad = pi/180.0
	cdef double rad_to_deg = 180.0/pi

	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180

	cdef double inv_LE = 1.0/L_E
	cdef double inv_LA = 1.0/L_A

	for i in range(N):
		lat1 = unkn_lat_arr[i]*deg_to_rad
		lon1 = unkn_lon_arr[i]*deg_to_rad
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for j in range(M):
			d_E  = vincenty_cy(data_lat_arr[j],data_lon_arr[j],unkn_lat_arr[i],unkn_lon_arr[i])
			rho_E       = exp(-pow(d_E*inv_LE,gamma_E))
			lat2        = data_lat_arr[j]*deg_to_rad
			lon2        = data_lon_arr[j]*deg_to_rad
			Az2         = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A         = acos(cos(Az1-Az2))*rad_to_deg
			rho_A       = (1+d_A*inv_LA)*pow(1-d_A/180,180*inv_LA)
			rho_EA      = rho_E*rho_A
			MrhoEA[j,i] = rho_EA
	return MrhoEA

#===================================================================================================
# rho_EAS: Sampled location to unsampled location
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEAS_sam2uns(data_lat, data_lon, data_Vs30, double epi_lat, double epi_lon, unkn_lat, unkn_lon, unkn_Vs30, double L_E, double gamma_E, double L_A, double L_S, double weight):

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

	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr  = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr  = np.ascontiguousarray(data_lon, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_Vs30_arr = np.ascontiguousarray(data_Vs30, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] unkn_lat_arr  = np.ascontiguousarray(unkn_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] unkn_lon_arr  = np.ascontiguousarray(unkn_lon, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] unkn_Vs30_arr = np.ascontiguousarray(unkn_Vs30, dtype=np.float64)

	cdef Py_ssize_t M = data_lat_arr.shape[0]
	cdef Py_ssize_t N = unkn_lat_arr.shape[0]

	cdef np.ndarray[np.double_t, ndim=2] MrhoEAS = np.empty((M, N), dtype=np.float64)

	cdef Py_ssize_t i, j
	cdef double lat1, lon1, Az1, d_E, rho_E, lat2, lon2, Az2, d_A, rho_A, d_S, rho_S, rho_EAS

	cdef double pi = 3.141592653589793
	cdef double deg_to_rad = pi/180.0
	cdef double rad_to_deg = 180.0/pi

	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180

	cdef double inv_LE = 1.0/L_E
	cdef double inv_LA = 1.0/L_A
	cdef double inv_LS = 1.0/L_S
	cdef double one_minus_w = 1-weight

	for i in range(N):
		lat1 = unkn_lat_arr[i]*deg_to_rad
		lon1 = unkn_lon_arr[i]*deg_to_rad
		Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
		for j in range(M):
			d_E          = vincenty_cy(data_lat_arr[j],data_lon_arr[j],unkn_lat_arr[i],unkn_lon_arr[i])
			rho_E        = exp(-pow(d_E*inv_LE,gamma_E))
			lat2         = data_lat_arr[j]*deg_to_rad
			lon2         = data_lon_arr[j]*deg_to_rad
			Az2          = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
			d_A          = acos(cos(Az1-Az2))*rad_to_deg
			rho_A        = (1+d_A*inv_LA)*pow(1-d_A/180,180*inv_LA)
			d_S          = fabs(unkn_Vs30_arr[i]-data_Vs30_arr[j])
			rho_S        = exp(-(d_S*inv_LS))
			rho_EAS      = rho_E*(weight*rho_A+one_minus_w*rho_S)
			MrhoEAS[j,i] = rho_EAS
	return MrhoEAS

#===================================================================================================
# rho_E: Sampled to unsampled grid
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoE_sam2grid(data_lat, data_lon, grid_lat, grid_lon, double L_E, double gamma_E):

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

	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr = np.ascontiguousarray(data_lon, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] grid_lat_arr = np.ascontiguousarray(grid_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] grid_lon_arr = np.ascontiguousarray(grid_lon, dtype=np.float64)

	cdef Py_ssize_t M      = data_lat_arr.shape[0]
	cdef Py_ssize_t n_lats = grid_lat_arr.shape[0]
	cdef Py_ssize_t n_lons = grid_lon_arr.shape[0]

	cdef np.ndarray[np.double_t, ndim=2] MrhoE = np.empty((M, n_lats*n_lons), dtype=np.float64)

	cdef Py_ssize_t i, j, k
	cdef double d_E, rho_E

	cdef double inv_LE = 1.0/L_E

	for i in range(n_lons):
		for j in range(n_lats):
			for k in range(M):
				d_E                  = vincenty_cy(data_lat_arr[k],data_lon_arr[k],grid_lat_arr[j],grid_lon_arr[i])
				rho_E                = exp(-pow(d_E*inv_LE,gamma_E))
				MrhoE[k,i*n_lats+j] = rho_E
	return MrhoE

#===================================================================================================
# rho_EA: Sampled to unsampled grid
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def MrhoEA_sam2grid(data_lat, data_lon, double epi_lat, double epi_lon, grid_lat, grid_lon, double L_E, double gamma_E, double L_A):

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


	cdef np.ndarray[np.double_t, ndim=1] data_lat_arr = np.ascontiguousarray(data_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] data_lon_arr = np.ascontiguousarray(data_lon, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] grid_lat_arr = np.ascontiguousarray(grid_lat, dtype=np.float64)
	cdef np.ndarray[np.double_t, ndim=1] grid_lon_arr = np.ascontiguousarray(grid_lon, dtype=np.float64)

	cdef Py_ssize_t M      = data_lat_arr.shape[0]
	cdef Py_ssize_t n_lats = grid_lat_arr.shape[0]
	cdef Py_ssize_t n_lons = grid_lon_arr.shape[0]

	cdef np.ndarray[np.double_t, ndim=2] MrhoEA = np.empty((M, n_lats*n_lons), dtype=np.float64)

	cdef Py_ssize_t i, j, k
	cdef double lat1, lat2, lon1, lon2, Az1, Az2, d_E, d_A, rho_E, rho_A, rho_EA

	cdef double pi = 3.141592653589793
	cdef double deg_to_rad = pi/180.0
	cdef double rad_to_deg = 180.0/pi

	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180

	cdef double inv_LE = 1.0/L_E
	cdef double inv_LA = 1.0/L_A

	for i in range(n_lons):
		for j in range(n_lats):
			lat1 = grid_lat_arr[j]*deg_to_rad
			lon1 = grid_lon_arr[i]*deg_to_rad
			Az1  = atan2(sin(lon1-epi_lon)*cos(lat1),cos(epi_lat)*sin(lat1)-sin(epi_lat)*cos(lat1)*cos(lon1-epi_lon))
			for k in range(M):
				d_E                  = vincenty_cy(data_lat_arr[k],data_lon_arr[k],grid_lat_arr[j],grid_lon_arr[i])
				rho_E                = exp(-pow(d_E*inv_LE,gamma_E))
				lat2                 = data_lat_arr[k]*deg_to_rad
				lon2                 = data_lon_arr[k]*deg_to_rad
				Az2                  = atan2(sin(lon2-epi_lon)*cos(lat2),cos(epi_lat)*sin(lat2)-sin(epi_lat)*cos(lat2)*cos(lon2-epi_lon))
				d_A                  = acos(cos(Az1-Az2))*rad_to_deg
				rho_A                = (1+d_A*inv_LA)*pow(1-d_A/180,180*inv_LA)
				rho_EA               = rho_E*rho_A
				MrhoEA[k,i*n_lats+j] = rho_EA
	return MrhoEA
