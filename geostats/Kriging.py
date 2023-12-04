""" Kriging tools.
"""
import numpy as np
import torch
from .build_Mrho import build_Mrho

__author__ = 'A. Renmin Pretell Ductram'


class Model:
	
	"""
	model_type: Correlation model, options are: 'E', 'EA', and 'EAS' (Bodenmann et al., 2023)
	L_E: Euclidean distance correlation length.
	gamma_E: Exponential coefficient in Euclidean correlation model.
	L_A: Azimuthal distance correlation length in deg.
	L_S: Vs30 dissimilarity distance correlation length.
	weight: Weight model coefficient.
	"""
	
	def __init__(self, model_type, L_E, gamma_E, L_A=None, L_S=None, weight=None):
		self.L_E     = L_E
		self.gamma_E = gamma_E
		
		if model_type.lower() == 'e':
			self.m_type = 'e'
		elif model_type.lower() == 'ea':
			self.L_A    = L_A
			self.m_type = 'ea'
		elif model_type.lower() == 'eas':
			self.L_A    = L_A
			self.L_S    = L_S
			self.weight = weight
			self.m_type = 'eas'


class Site:
	
	"""
	latitude: Sites latitude in deg.
	longitude: Sites longitude in deg.
	Vs30: Time-averaged shear wave velocity in the top 30 m.
	value: Values for parameter to be Kriged.
	"""
	
	def __init__(self, latitude, longitude, Vs30=None, value=None):
		self.value     = value
		self.latitude  = latitude
		self.longitude = longitude
		self.Vs30      = Vs30


class Kriging:
	
	"""
	sam: Data points, normalized to have a variance of 1. Must be a 'site' instance.
	uns: Unknown points. Must be a 'site' instance.
	var: Data variance.
	model: Correlation model. Must be a 'model' instance.
	epi_lat/epi_lon: Coordinates of the epicenter in deg.
	Mrho_sam: Correlation matrix for data points. 
	Mrho_uns: Correlation matrix between data and unknown points. 
	"""
	
	def __init__(self, sam, uns, var, model=None, epi_lat=None, epi_lon=None, Mrho_sam=None, Mrho_uns=None):
		self.sam   = sam
		self.uns   = uns
		self.model = model
		self.var   = var
		
		if epi_lat is not None:
			self.epi_lat = epi_lat
			self.epi_lon = epi_lon
		
		if Mrho_sam is not None:
			self.Mrho_sam = Mrho_sam
			self.Mrho_uns = Mrho_uns
		else:
			self.Mrho_sam = None
			self.Mrho_uns = None
	
	class Kriged:
		
		def __init__(self,mean,variance):
			self.mean     = mean
			self.variance = variance

	def Krige(self,Ktype):
		
		if self.Mrho_sam is None:
		
			if self.model.m_type.lower() == 'e':
				self.Mrho_sam = build_Mrho.MrhoE_sam2sam_cy(self.sam.latitude,self.sam.longitude,self.model.L_E,self.model.gamma_E)
				self.Mrho_uns = build_Mrho.MrhoE_sam2uns_cy(self.sam.latitude,self.sam.longitude,self.uns.latitude,self.uns.longitude,self.model.L_E,self.model.gamma_E)
			
			elif self.model.m_type.lower() == 'ea':
				self.Mrho_sam = build_Mrho.MrhoEA_sam2sam_cy(self.sam.latitude,self.sam.longitude,self.epi_lat,self.epi_lon,self.model.L_E,self.model.gamma_E,self.model.L_A)
				self.Mrho_uns = build_Mrho.MrhoEA_sam2uns_cy(self.sam.latitude,self.sam.longitude,self.epi_lat,self.epi_lon,self.uns.latitude,self.uns.longitude,self.model.L_E,self.model.gamma_E,self.model.L_A)
		
			elif self.model.m_type.lower() == 'eas':
				self.Mrho_sam = build_Mrho.MrhoEAS_sam2sam_cy(self.sam.latitude,self.sam.longitude,self.sam.Vs30,self.epi_lat,self.epi_lon,self.model.L_E,self.model.gamma_E,self.model.L_A,self.model.L_S,self.model.weight)
				self.Mrho_uns = build_Mrho.MrhoEAS_sam2uns_cy(self.sam.latitude,self.sam.longitude,self.sam.Vs30,self.epi_lat,self.epi_lon,self.uns.latitude,self.uns.longitude,self.uns.Vs30,self.model.L_E,self.model.gamma_E,self.model.L_A,self.model.L_S,self.model.weight)
	
		M = len(self.Mrho_sam)
		N = len(self.Mrho_uns.T)
		
		Mcov_sam_ = self.Mrho_sam*self.var
		Mcov_uns  = self.Mrho_uns*self.var
		
		if Ktype.lower() == 'ordinary':
			Mcov_uns = np.vstack((Mcov_uns,np.ones(N)))
			Mcov_sam = np.hstack((Mcov_sam_,np.transpose([np.ones(M)])))
			Mcov_sam = np.vstack((Mcov_sam,np.ones(M+1)))
			Mcov_sam[M,M] = 0.0
			A = torch.tensor(Mcov_sam)
			b = torch.tensor(Mcov_uns)
			k = torch.linalg.solve(A,b)
			w  = k[0:M,0:N].T
			mu = w@self.sam.value
			vr = self.var-(w@Mcov_uns[0:M])-k.T[:,-1]
			vr = np.diagonal(vr)
		else:
			raise NotImplementedError
		
		return self.Kriged(mu.numpy(),vr)


