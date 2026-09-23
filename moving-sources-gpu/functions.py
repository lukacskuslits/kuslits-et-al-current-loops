import os
import pandas as pd
import csv
import time
from datetime import datetime

import torch
import shtns
import numpy as np

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


# ============================================================
# Jacobi polynomials (vectorized, no tensor recreation)
# ============================================================
def jacobi_recurrence_for_GPU(power, alpha, beta, x):
    if power < 0:
        return torch.zeros_like(x)

    J0 = torch.ones_like(x)
    if power == 0:
        return J0

    J1 = 0.5 * (alpha + beta + 2) * x + 0.5 * (alpha - beta)
    if power == 1:
        return J1

    J_prev = J0
    J_curr = J1

    for n in range(2, power + 1):
        an = (2*(n-1)+alpha+beta+1)*(2*(n-1)+alpha+beta+2) / (2*n*(n+alpha+beta))
        bn = (beta**2 - alpha**2)*(2*(n-1)+alpha+beta+1) / (2*n*(n+alpha+beta)*(2*(n-1)+alpha+beta))
        cn = (n+1+alpha)*(n+1+beta)*(2*(n+1)+alpha+beta+2) / ((n+2)*(n+alpha+beta+2)*(2*(n+1)+alpha+beta))

        J_next = (an * x - bn) * J_curr - cn * J_prev
        J_prev = J_curr
        J_curr = J_next

    return J_curr


# ============================================================
# Psi computation (vectorized)
# ============================================================
def compute_all_Psi(r, degree, band_nr, mu):
    r = r.to(device)
    z = 2 * r**2 - 1

    alpha = 2
    beta = degree + 0.5

    Psi_vals = []

    for k in range(1, band_nr):
        if k == 1:
            Psi = r**degree * (-mu*degree - 2*mu - 1 + r**2 * mu*degree + r**2)
        elif k == 2:
            Psi = r**degree * (
                6*mu*degree**2 + 27*mu*degree + 27*mu + 2*degree + 3
                - 12*r**2*degree**2*mu - 4*r**2*degree
                - 58*r**2*mu*degree - 10*r**2 - 70*r**2*mu
                + 6*r**4*degree**2*mu + 2*r**4*degree
                + 31*r**4*mu*degree + 35*r**4*mu + 7*r**4
            )
        else:
            c = torch.zeros(4, device=device)

            c[0] = (2*degree - 3 + 4*k)*(2*degree + 5 + 2*k)*(2*degree - 1 + 4*k)*(2*k*mu*degree - mu*degree + 1 + 2*k**2*mu - mu - k*mu)*(2*degree + 3 + 2*k)
            c[1] = -(2*degree - 3 + 4*k)*(2*degree + 3 + 2*k)*(2*degree + 3 + 4*k)*(2*degree + 1 + 2*k)*(-mu*degree + 6*k*mu*degree + 3 + k*mu + 6*k**2*mu + 2*mu)
            c[2] = (2*degree - 1 + 4*k)*(2*degree + 5 + 4*k)*(2*degree + 1 + 2*k)*(-1 + 2*degree + 2*k)*(6*k*mu*degree + mu*degree + 5*k*mu + 3 + 3*mu + 6*k**2*mu)
            c[3] = -(2*degree - 3 + 2*k)*(2*degree + 5 + 4*k)*(2*degree + 3 + 4*k)*(-1 + 2*degree + 2*k)*(2*k*mu*degree + mu*degree + 3*k*mu + 2*k**2*mu + 1)

            jacobi_vec = torch.stack([
                jacobi_recurrence_for_GPU(k, alpha, beta, z),
                jacobi_recurrence_for_GPU(k-1, alpha, beta, z),
                jacobi_recurrence_for_GPU(k-2, alpha, beta, z),
                jacobi_recurrence_for_GPU(k-3, alpha, beta, z)
            ])

            Psi = r**degree * torch.matmul(c, jacobi_vec)

        Psi_vals.append(Psi)

    return torch.stack(Psi_vals)  # shape: [band_nr-1, N]


# ============================================================
# Vectorized integration (NO loops, NO torchquad)
# ============================================================
def discretized_laplacian_B_fast(degree, r0, r1, resolution, band_nr, normalization_mx, l2):
    mu = 1 / l2[degree]

    r = torch.linspace(r0, r1, resolution, device=device)
    dr = (r1 - r0) / (resolution - 1)

    Psi_vals = compute_all_Psi(r, degree, band_nr, mu)

    # normalization
    norm = normalization_mx[:, degree].to(device)
    Psi_vals = Psi_vals * norm[:, None]

    weight = r**2

    # Compute all integrals at once
    B = torch.einsum('ir,jr,r->ij', Psi_vals, Psi_vals, weight) * dr

    return B


# ============================================================
# Normalization constants (slightly optimized)
# ============================================================
def derive_normalization_constant(r, band_nr, max_degree, l2):
    r = r.to(device)
    NS = torch.ones([band_nr-1, max_degree], device=device)

    for degree in range(1, max_degree):
        mu = 1 / l2[degree]
        Psi_vals = compute_all_Psi(r.squeeze(), degree, band_nr, mu)
        NS[:, degree] = 1/torch.max(torch.abs(Psi_vals), dim=1).values

    NS[NS == 0] = 1
    return NS


# ============================================================
# Jacobi polynomials (vectorized, no tensor recreation)
# ============================================================
def jacobi_recurrence_for_GPU(power, alpha, beta, x):
    if power < 0:
        return torch.zeros_like(x)

    J0 = torch.ones_like(x)
    if power == 0:
        return J0

    J1 = 0.5 * (alpha + beta + 2) * x + 0.5 * (alpha - beta)
    if power == 1:
        return J1

    J_prev = J0
    J_curr = J1

    for n in range(2, power + 1):
        an = (2*(n-1)+alpha+beta+1)*(2*(n-1)+alpha+beta+2) / (2*n*(n+alpha+beta))
        bn = (beta**2 - alpha**2)*(2*(n-1)+alpha+beta+1) / (2*n*(n+alpha+beta)*(2*(n-1)+alpha+beta))
        cn = (n+1+alpha)*(n+1+beta)*(2*(n+1)+alpha+beta+2) / ((n+2)*(n+alpha+beta+2)*(2*(n+1)+alpha+beta))

        J_next = (an * x - bn) * J_curr - cn * J_prev
        J_prev = J_curr
        J_curr = J_next

    return J_curr


# ============================================================
# Psi computation (vectorized)
# ============================================================
def compute_all_Psi(r, degree, band_nr, mu):
    r = r.to(device)
    z = 2 * r**2 - 1

    alpha = 2
    beta = degree + 0.5

    Psi_vals = []

    for k in range(1, band_nr):
        if k == 1:
            Psi = r**degree * (-mu*degree - 2*mu - 1 + r**2 * mu*degree + r**2)
        elif k == 2:
            Psi = r**degree * (
                6*mu*degree**2 + 27*mu*degree + 27*mu + 2*degree + 3
                - 12*r**2*degree**2*mu - 4*r**2*degree
                - 58*r**2*mu*degree - 10*r**2 - 70*r**2*mu
                + 6*r**4*degree**2*mu + 2*r**4*degree
                + 31*r**4*mu*degree + 35*r**4*mu + 7*r**4
            )
        else:
            c = torch.zeros(4, device=device)

            c[0] = (2*degree - 3 + 4*k)*(2*degree + 5 + 2*k)*(2*degree - 1 + 4*k)*(2*k*mu*degree - mu*degree + 1 + 2*k**2*mu - mu - k*mu)*(2*degree + 3 + 2*k)
            c[1] = -(2*degree - 3 + 4*k)*(2*degree + 3 + 2*k)*(2*degree + 3 + 4*k)*(2*degree + 1 + 2*k)*(-mu*degree + 6*k*mu*degree + 3 + k*mu + 6*k**2*mu + 2*mu)
            c[2] = (2*degree - 1 + 4*k)*(2*degree + 5 + 4*k)*(2*degree + 1 + 2*k)*(-1 + 2*degree + 2*k)*(6*k*mu*degree + mu*degree + 5*k*mu + 3 + 3*mu + 6*k**2*mu)
            c[3] = -(2*degree - 3 + 2*k)*(2*degree + 5 + 4*k)*(2*degree + 3 + 4*k)*(-1 + 2*degree + 2*k)*(2*k*mu*degree + mu*degree + 3*k*mu + 2*k**2*mu + 1)

            jacobi_vec = torch.stack([
                jacobi_recurrence_for_GPU(k, alpha, beta, z),
                jacobi_recurrence_for_GPU(k-1, alpha, beta, z),
                jacobi_recurrence_for_GPU(k-2, alpha, beta, z),
                jacobi_recurrence_for_GPU(k-3, alpha, beta, z)
            ])

            Psi = r**degree * torch.matmul(c, jacobi_vec)

        Psi_vals.append(Psi)

    return torch.stack(Psi_vals)  # shape: [band_nr-1, N]


# ============================================================
# Vectorized integration (NO loops, NO torchquad)
# ============================================================
def discretized_laplacian_B_fast(degree, r0, r1, resolution, band_nr, normalization_mx, l2):
    mu = 1 / l2[degree]

    r = torch.linspace(r0, r1, resolution, device=device)
    dr = (r1 - r0) / (resolution - 1)

    Psi_vals = compute_all_Psi(r, degree, band_nr, mu)

    # normalization
    norm = normalization_mx[:, degree].to(device)
    Psi_vals = Psi_vals * norm[:, None]

    weight = r**2

    # Compute all integrals at once
    B = torch.einsum('ir,jr,r->ij', Psi_vals, Psi_vals, weight) * dr

    return B


# ============================================================
# Normalization constants (slightly optimized)
# ============================================================
def derive_normalization_constant(r, band_nr, max_degree, l2):
    r = r.to(device)
    NS = torch.ones([band_nr-1, max_degree], device=device)

    for degree in range(1, max_degree):
        mu = 1 / l2[degree]
        Psi_vals = compute_all_Psi(r.squeeze(), degree, band_nr, mu)
        NS[:, degree] = 1/torch.max(torch.abs(Psi_vals), dim=1).values

    NS[NS == 0] = 1
    return NS

def compute_toroidal_flow_coefficients_vectorized(n, m, r, Omega):
    phi = torch.Tensor([torch.pi / 2]).to(device)
    th = torch.Tensor([torch.pi / 4]).to(device)
    #for th in list(torch.linspace(0.01, torch.pi/2-0.01, 1000, device=device)):
    #TODO: DOES THE VALUE OF T_10 REMAIN THE SAME FOR EVERY PHI?
    t_10 = -Omega*torch.sin(th)*r/Y_derivative_in_theta(n, m, th, phi)
    t_10 = t_10[~torch.isnan(t_10)]
    #TODO DOUBLE CHECK THIS TO STABILIZE VALUE:
    return t_10

#-------------------------- COMPUTE THE COEFFICIENT OF A PRESCRIBED TOROIDAL FLOW --------------------------------------
def compute_toroidal_flow_coefficient(n, m, r, Omega):
    t_10 = torch.empty([0,]).to(device)
    phi = torch.pi / 2
    th = torch.pi / 4
    for th in list(torch.linspace(0.01, torch.pi/2-0.01, 1000, device=device)):
    #TODO: DOES THE VALUE OF T_10 REMAIN THE SAME FOR EVERY PHI?
        t_10 = torch.cat((t_10, -Omega*r*torch.sin(th)/Y_derivative_in_theta(n, m, th, phi)),0)
    t_10 = t_10[~torch.isnan(t_10)]
    #TODO DOUBLE CHECK THIS TO STABILIZE VALUE:
    t_10 = t_10[500]
    return t_10

def define_3d_grid_for_vector_SH_transform(sh):
    x = sh.spat_array()  # a spatial array, same as numpy.zeros(sh.spat_shape)
    y = x.copy()
    y[:, :] = 0
    z = x.copy()
    z[:, :] = 0
    return x, y, z


def compute_poloidal_source_coefficients(Br_test, x, y, z, sh, r_samples, r_ind):
    """
    Performs the spherical harmonic analysis on a map of source field values using an instance of the shtns package (sh)
    with three arguments, sh.alanys() performs a 3D vector transform (spherical coordinates)
    :param Br_test: stack of geographic maps (radial slices) of source field values
    :param x: starting grid
    :param y:
    :param z:
    :param sh: an instance of the shtns package
    :param r_samples: radial sampling points for geomagnetic map slices
    :param r_ind: index of the given geocentric radius of a slice
    :return: S0_r source coefficients; el array of size sh.nlm giving the SH degree l for any SH coefficient,
    l2 array l(l+1) that is useful for computing laplacian
    """
    x[:, :] = torch.squeeze(Br_test[:, :, r_ind]).cpu().numpy()
    qlm, slm, tlm = sh.analys(x, y,
                              z)
    el = sh.l
    em = sh.m
    l2 = el * (el + 1)
    r_dimensional = r_samples[r_ind]
    qlm = torch.Tensor(qlm).to(device)
    l2 = torch.Tensor(l2).to(device)
    S0_r = (qlm * r_dimensional) / (l2 + 1)
    return S0_r, el, em, l2

def create_vector_of_galerkin_polynomials(mu, gamma_two_index, normalization_mx, band_nr, r, el):
    """
    Constructing the Galerkin-base at each observation radius containing basis functions of each polynomial degree used
     for inverting the coefficients of the Galerkin-representation of the induced field
    :param mu: degree-based multiplier
    :param gamma_two_index: SH two-index gamma of the induced field
    :param normalization_mx: array A containing the normalization constants
    :param band_nr: maximum polynomial degree used for the representation
    :param r: radial observation point
    :param el: array containing all SH degrees
    :return: [A^1_{\gamma}Psi^1_{\gamma}(r), A^2_{\gamma}Psi^2_{\gamma}(r), ..., A^K_{\gamma}Psi^K_{\gamma}(r)]
    """
    psi_val = compute_all_Psi(r, el[gamma_two_index], band_nr, mu)
    psi_vec = torch.mul(torch.squeeze(psi_val), torch.squeeze(normalization_mx[:,gamma_two_index].T))
    # psi_vec = torch.empty([0,]).to(device)
    # for jacobi_k in range(band_nr):
    #     psi_val = compute_all_Psi(r, el[gamma_two_index], band_nr, mu)[jacobi_k]
    #     psi_vec = torch.cat((psi_vec, (1 / normalization_mx[jacobi_k, gamma_two_index]) * psi_val), 0)
    #psi_vec = psi_vec[1:]
    return psi_vec

# -------------------------ELSASSER DYNAMO INTEGRAL AS FEATURED IN [5] --------------------------------------------------
def elsasser(n,m,k,l,i,j):
    """
    Evaluates the Elsasser dynamo integral for each corresponding SH degree and order listed below
    :param n: flow order n_gamma
    :param m: flow degree m_gamma
    :param k: source mf order k_beta
    :param l: source mf degree l_beta
    :param i: induced mf order i_gamma
    :param j: induced mf degree j_gamma
    :return: elsasser_nm value of the Elsasser dynamo integral
    """
    Lambda_b_g = cmath.sqrt((2 * n + 1)*(2 * k + 1)*(2 * i + 1))
    Delta_b_g = cmath.sqrt((n + k + i + 2)*(n + k + i + 4) / (4 * (n + k + i + 3))) * cmath.sqrt((n + k - i + 1)*(n + i - k + 1)*(i + k - n + 1))
    elsasser_nm = (-1)**np.float32(j)*(-4j*pi*Lambda_b_g*Delta_b_g*wigner_3j(n+1, k+1, i + 1, 0, 0, 0)*wigner_3j(n, k, i, m, l, j))
    return elsasser_nm

