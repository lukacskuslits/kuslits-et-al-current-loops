import numpy as np

n_gamma = 11
p = 0.5
rho = 1.25
r = 3.48e6

Omega = 5*1.2e-10
mu_0 = 4*np.pi*1e-7
sigma = 1e6
eta = 1/(mu_0*sigma)
#TODO: vary elsasser_gamma with n_gamma!!
elsasser_gamma = np.array(21.77, dtype='complex')


lambda_gamma_roots = {'n_gamma': [7, 8, 9, 10, 11], 'lambda_gamma': [66, 72, 82, 90, 98]}