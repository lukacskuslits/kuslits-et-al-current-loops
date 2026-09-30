import numpy as np
from pyshtools.legendre import PlmScmidt

def alldredge_loop_sh_coeffs(n, m, K, alpha, r, theta_0, phi_0):
	'''
	# Computes the sine and cosine SHCs of order n, degree m for a single current loop
	# Input parameters:
	#-----------------------
	n: SH degree
	m: SH order
	K: normalized magnetic moment
	alpha: half angle of sight from EC to the loop
	r: distance from EC to the edge of the loop
	theta_0: spherical colatitude
	phi_0: spherical longitude
	Output values:
	----------------------
	gnm sine SHC
	hnm cosine SHC
	'''

	a = 6.378e6 # Earth's radius [m]

	x = np.cos(alpha)
	legendre_1 = PlmScmidt(n, x, [1,0])[0]
	g0n = K*((2*n)/(np.sin(alpha)))*((r/a)^(n-1))*(legendre_1/np.sqrt(2*n*(n+1)))

	x_m = np.cos(theta_0)
	legendre_m = PlmScmidt(n, x_m, [1,0])[0]
	gnm = g0n*legendre_m*np.cos(m*phi_0)
	hnm = g0n*legendre_m*np.sin(m*phi_0)
	return gnm, hnm