import torch
import numpy as np
import shtns
import scipy

from functions import *


# ============================================================
# GLOBAL_DEFINITIONS
# ============================================================

global mu_0
mu_0 = 4 * np.pi * 1e-7
sigma = 5e5
global eta
eta = 1 / (mu_0 * sigma)

lmax = 14
mmax = 14

lat_max = 90
long_max = 180

global sh
sh = shtns.sht(lmax, mmax)
nlat, nphi = sh.set_grid(lat_max, long_max)

global el
el = sh.l
global em
em = sh.m
global l2
l2 = el * (el + 1)

global band_nr
band_nr = 13
degree = 14
global r0
r0 = 0.0
global r1
r1 = 1.0
global resolution
resolution = 1000

# Number of radial sample (observation) points
global nr_r_samples
nr_r_samples = 25

# Sample radii
global r_samples
r_samples = torch.linspace(0.01, 1, nr_r_samples, device=device)

# Maximum primary source mf two-indexes
global beta_max
beta_max = 120

# Maximum induced mf two-indexes
global gamma_max
gamma_max = 120

# ============================================================
# NORMALIZATION
# ============================================================
global normalization_mx
normalization_mx = derive_normalization_constant(r_samples, band_nr, gamma_max, l2)

# Angular velocity of rotating flow
global Omega
Omega = 0.01  # [1/s]
# Observation radius (normalized to the CMB)
r = 0.95
# Alpha indices - for the toroidal flow coefficient - FIXED
global n
n = 1
global m
m = 0
# Matrix containing the full Galerkin base $Psi(r)$
psi_tmp_local = torch.ones([25, band_nr-1, gamma_max]).to(device)
r_ind = 0
for r in r_samples:
    for gamma_two_index in range(1, gamma_max):
        mu = 1 / l2[gamma_two_index]
        psi_tmp_local[r_ind, :, gamma_two_index] = create_vector_of_galerkin_polynomials(mu, gamma_two_index, normalization_mx, band_nr, r, el)
    r_ind += 1
global PSI_MATRIX
PSI_MATRIX = psi_tmp_local
# Array containing magnetic field values
global Br_test
Br_test = torch.Tensor(np.array(scipy.io.loadmat('B_r_test.mat')['B_r'])).to(device)
# First SH two index of the GB solution iteration
global gamma_first_index
gamma_first_index = 37

def forward_computation_BG_residual(q):
    print("Forward function called with:", flush=True)
    print(q, flush=True)
    total_res = torch.Tensor([0]).to(device)
    for gamma_two_index in range(gamma_first_index, globals()['gamma_max']):
        ## COMPUTE THE RHS OF THE B-G EQs ##########################################################
        RHS_r = torch.empty([0, ]).to(device)
        LHS_r = torch.empty([0, ]).to(device)
        i = globals()['el'][gamma_two_index]
        j = globals()['em'][gamma_two_index]
        mu = 1 / (i + 1)
        t_10 = compute_toroidal_flow_coefficients_vectorized(n, m, globals()['r_samples'], globals()['Omega'])
        r_ind = 0
        S_induced = torch.zeros([25, gamma_max]).to(device)
        for r in r_samples:
            nu_rhs = 0
            x, y, z = define_3d_grid_for_vector_SH_transform(globals()['sh'])
            S0_r, el, em, l2 = compute_poloidal_source_coefficients(globals()['Br_test'], x, y, z, globals()['sh'],
                                                                    globals()['r_samples'], r_ind)
            for beta_two_index in list(range(1, globals()['beta_max'])):
                k = globals()['el'][beta_two_index]
                l = globals()['em'][beta_two_index]
                # psi_vec = create_vector_of_galerkin_polynomials(mu, beta_two_index, normalization_mx, band_nr, r, el)
                psi_vec = torch.squeeze(globals()['PSI_MATRIX'][r_ind, :, beta_two_index])
                S_induced[r_ind, beta_two_index] = torch.matmul(psi_vec, q.T)
                S_kl_r = S0_r[beta_two_index] + S_induced[r_ind, beta_two_index]
                Elsasser = elsasser(n, m, k, l, i, -j)
                Elsasser = np.abs(Elsasser)
                nu_rhs += (k * (k + 1)) * t_10[r_ind] * S_kl_r * Elsasser / (4 * np.pi * i * (i + 1))
            # TODO: RE-CHECK WHETHER THE LENGTH OR THE REAL PART OF THE COMPLEX VALUES NEED TO BE USED HERE:
            nu_rhs = torch.Tensor([nu_rhs]).to(device)
            RHS_r = torch.cat((RHS_r, nu_rhs), 0)
            r_ind += 1

        # COMPUTE THE LHS OF THE BG EQs ########################################
        B = discretized_laplacian_B_fast(
            globals()['el'][gamma_two_index],
            globals()['r0'],
            globals()['r1'],
            globals()['resolution'],
            globals()['band_nr'],
            globals()['normalization_mx'],
            l2
        )
        ##################
        E_matrix = torch.zeros([nr_r_samples, globals()['band_nr'] - 1]).to(device)
        ##################
        r_ind = 0
        for r in globals()['r_samples']:
            psi_vec = torch.squeeze(globals()['PSI_MATRIX'][r_ind, :, gamma_two_index])
            E_matrix[r_ind, :] = -eta * torch.matmul(torch.linalg.inv(B), psi_vec).T
            r_ind += 1
        LHS_r = torch.matmul(E_matrix, q)

        residual = torch.abs(LHS_r - RHS_r)
        total_res += torch.sum(residual)
    return total_res.cpu().numpy()


