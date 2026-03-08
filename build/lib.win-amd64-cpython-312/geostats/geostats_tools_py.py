""" A suite of functions supporting geotatistical analysis.
"""

import numpy as np
from vincenty import vincenty

__author__ = 'A. Renmin Pretell Ductram'

#===================================================================================================
# Euclidean epicentral distance
#===================================================================================================
def get_E_epi_dist(lats, lons, epi_lat, epi_lon):

	"""
	Parameter
	=========
	lats/lons: Site coordinates in deg.
	epi_lat/epi_lon: Epicenter coordinates in deg.
	
	Returns
	=======
	Euclidean epicentral distance in km.
	"""

	lats_arr = np.atleast_1d(lats)
	lons_arr = np.atleast_1d(lons)

	n_sites    = lats_arr.size
	epi_E_dist = np.empty(n_sites, dtype=np.float64)

	for i in range(n_sites):
		epi_E_dist[i] = vincenty((epi_lat,epi_lon),(lats_arr[i],lons_arr[i]))
	return epi_E_dist

#===================================================================================================
# Azimuthal distance
#===================================================================================================
# Modified after Lukas Bodenmann: https://github.com/bodlukas/ground-motion-correlation-bayes
def get_A_epi_dist(lats, lons, epi_lat, epi_lon):

	"""
	Parameter
	=========
	lats/lons: Site coordinates in deg.
	epi_lat/epi_lon: Coordinates of the epicenter in deg.
	
	Returns
	=======
	Azimuthal epicentral distance in rad.
	"""

	lats_arr = np.atleast_1d(lats)
	lons_arr = np.atleast_1d(lons)

	lats_arr_r = np.deg2rad(lats_arr)
	lons_arr_r = np.deg2rad(lons_arr)

	epi_lat = np.deg2rad(epi_lat)
	epi_lon = np.deg2rad(epi_lon)

	return np.atan2(np.sin(lons_arr_r-epi_lon)*np.cos(lats_arr_r),np.cos(epi_lat)*np.sin(lats_arr_r)-np.sin(epi_lat)*np.cos(lats_arr_r)*np.cos(lons_arr_r-epi_lon))

#===================================================================================================
# Vs30 dissimilarity
#===================================================================================================
def get_Vs30_dist(site_Vs30):

	"""
	Parameter
	=========
	site_Vs30: Site Vs30.
	
	Returns
	=======
	Vs30 dissimilarity values.
	"""

	site_Vs30_arr = np.asarray(site_Vs30, dtype=np.float64)
	return np.abs(site_Vs30_arr[:, None] - site_Vs30_arr[None, :])

#===================================================================================================
# Euclidean distance
#===================================================================================================
def get_E_dist(d_E, d_A):

	"""
	Parameters
	==========
	d_E: Epicentral Euclidean distances.
	d_A: Epicentral azimuthal distances in rad.
	
	Returns
	=======
	Euclidean distance matrix in km.
	"""

	d_E_arr = np.asarray(d_E, dtype=np.float64)
	d_A_arr = np.asarray(d_A, dtype=np.float64)

	dEi = d_E_arr[:, None]
	dEj = d_E_arr[None, :]

	dAi = d_A_arr[:, None]
	dAj = d_A_arr[None, :]

	return np.sqrt(dEi**2 + dEj**2 - 2.0 * dEi * dEj * np.cos(dAi - dAj))

#===================================================================================================
# Azimuthal distance
#===================================================================================================
def get_A_dist(d_A):

	"""
	Parameters
	=========
	d_A: Epicentral azimuthal distances in rad.
	
	Returns
	=======
	Azimuthal distance matrix in rad.
	"""

	d_A_arr = np.asarray(d_A, dtype=np.float64)

	delta_dA = d_A_arr[:, None] - d_A_arr[None, :]

	return np.arccos(np.cos(delta_dA))
