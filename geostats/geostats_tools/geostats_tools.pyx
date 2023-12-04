""" A suite of functions supporting geotatistical analysis.
"""

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport pi, sin, cos, tan, atan, atan2, acos, pow, sqrt

__author__ = 'A. Renmin Pretell Ductram'

#===================================================================================================
# Distance between two locations on the earth
#===================================================================================================
# Cythonized version of the vicenty function by Maurycy Pietrzak: https://github.com/maurycyp/vincenty
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def vincenty_cy(double lat1, double lon1, double lat2, double lon2):

	"""
	Parameter
	=========
	lat/lon: Site coordinates in deg.
	
	Returns
	=======
	Distance in km.
	"""
	
	cdef int    a = 6378137  # meters
	cdef double f = 1 / 298.257223563
	cdef double b = 6356752.314245
	cdef int max_iter = 200
	cdef double conv  = 1e-12
	cdef double sinLambda, cosLambda, sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, C, LambdaPrev

	if (lat1-lat2) == 0 and (lon1-lon2) == 0:
		return 0.0

	lat1 = lat1*pi/180
	lat2 = lat2*pi/180
	lon1 = lon1*pi/180
	lon2 = lon2*pi/180

	cdef double U1     = atan((1-f)*tan(lat1))
	cdef double U2     = atan((1-f)*tan(lat2))
	cdef double Lambda = lon2-lon1
	cdef double L      = lon2-lon1
	cdef double sinU1  = sin(U1)
	cdef double cosU1  = cos(U1)
	cdef double sinU2  = sin(U2)
	cdef double cosU2  = cos(U2)

	for iteration in range(max_iter):
		sinLambda = sin(Lambda)
		cosLambda = cos(Lambda)
		sinSigma  = sqrt(pow(cosU2*sinLambda,2)+pow(cosU1*sinU2-sinU1*cosU2*cosLambda,2))
		
		if sinSigma == 0:
			return 0.0
		
		cosSigma   = sinU1*sinU2+cosU1*cosU2*cosLambda
		sigma      = atan2(sinSigma,cosSigma)
		sinAlpha   = cosU1*cosU2*sinLambda/sinSigma
		cosSqAlpha = 1-pow(sinAlpha,2)
		
		try:
			cos2SigmaM = cosSigma-2*sinU1*sinU2/cosSqAlpha
		except ZeroDivisionError:
			cos2SigmaM = 0
		C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
		LambdaPrev = Lambda
		Lambda     = L+(1-C)*f*sinAlpha*(sigma+C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*pow(cos2SigmaM,2))))
		
		if abs(Lambda-LambdaPrev)<conv:
			break
	else:
		return None

	cdef double uSq    = cosSqAlpha*(pow(a,2)-pow(b,2))/pow(b,2)
	cdef double A      = 1+uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
	cdef double B      = uSq/1024*(256+uSq*(-128+uSq*(74-47*uSq)))
	cdef double dSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*pow(cos2SigmaM,2))-B/6*cos2SigmaM*(-3+4*pow(sinSigma,2))*(-3+4*pow(cos2SigmaM,2))))
	cdef double s      = b*A*(sigma-dSigma)/1000

	return s

#===================================================================================================
# Euclidean epicentral distance
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def get_E_epi_dist_cy(double[:] lats,double[:] lons,double epi_lat,double epi_lon):

	"""
	Parameter
	=========
	lats/lons: Site coordinates in deg.
	epi_lat/epi_lon: Epicenter coordinates in deg.
	
	Returns
	=======
	Euclidean epicentral distance in km.
	"""
	
	cdef int n_sites = int(len(lats))
	cdef int i
	cdef double[:] epi_E_dist = np.zeros([n_sites],dtype='float64')

	for i in range(n_sites):
		epi_E_dist[i] = vincenty_cy(epi_lat,epi_lon,lats[i],lons[i])
	return np.asarray(epi_E_dist)

#===================================================================================================
# Azimuthal distance
#===================================================================================================
# Modified after Lukas Bodenmann: https://github.com/bodlukas/ground-motion-correlation-bayes
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def get_A_epi_dist_cy(double[:] lats, double[:] lons, double epi_lat, double epi_lon):

	"""
	Parameter
	=========
	lats/lons: Site coordinates in deg.
	epi_lat/epi_lon: Coordinates of the epicenter in deg.
	
	Returns
	=======
	Azimuthal epicentral distance in rad.
	"""
	
	cdef int n_sites = int(len(lats))
	cdef int i
	cdef double[:] epi_A_dist = np.zeros([n_sites],dtype='float64')
	
	epi_lat = epi_lat*pi/180
	epi_lon = epi_lon*pi/180
	for i in range(n_sites):
		lat = lats[i]*pi/180
		lon = lons[i]*pi/180
		epi_A_dist[i] = atan2(sin(lon-epi_lon)*cos(lat),cos(epi_lat)*sin(lat)-sin(epi_lat)*cos(lat)*cos(lon-epi_lon))
	return np.asarray(epi_A_dist)

#===================================================================================================
# Vs30 dissimilarity
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def get_Vs30_dist_cy(double[:] site_Vs30):

	"""
	Parameter
	=========
	site_Vs30: Site Vs30.
	
	Returns
	=======
	Vs30 dissimilarity values.
	"""
	
	cdef int n_sites = int(len(site_Vs30))
	cdef int i, j
	cdef double[:,:] Vs30_dist = np.zeros([n_sites,n_sites],dtype='float64')

	for i in range(n_sites):
		for j in range(n_sites):
			if i != j:
				Vs30_dist[i][j] = abs(site_Vs30[i]-site_Vs30[j])
	return np.asarray(Vs30_dist)

#===================================================================================================
# Euclidean distance
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def get_E_dist_cy(double[:] d_E, double[:] d_A):
	
	"""
	Parameters
	==========
	d_E: Epicentral Euclidean distances.
	d_A: Epicentral azimuthal distances in rad.
	
	Returns
	=======
	Euclidean distance matrix in km.
	"""

	cdef int n_sites = int(len(d_E))
	cdef double[:,:] E_dist = np.zeros([n_sites,n_sites],dtype='float64')
	cdef int i, j

	for i in range(n_sites):
		for j in range(n_sites):
			if i != j:
				E_dist[i][j] = pow(d_E[i]**2+d_E[j]**2-2*d_E[i]*d_E[j]*cos(abs(d_A[i]-d_A[j])),0.5)
	return np.asarray(E_dist)

#===================================================================================================
# Azimuthal distance
#===================================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def get_A_dist_cy(double[:] d_A):
	
	"""
	Parameters
	=========
	d_A: Epicentral azimuthal distances in rad.
	
	Returns
	=======
	Azimuthal distance matrix in rad.
	"""
	
	cdef int n_sites = int(len(d_A))
	cdef double[:,:] A_dist = np.zeros([n_sites,n_sites],dtype='float64')
	cdef int i, j

	for i in range(n_sites):
		for j in range(n_sites):
			A_dist[i][j] = acos(cos(d_A[i]-d_A[j]))
	return np.asarray(A_dist)
