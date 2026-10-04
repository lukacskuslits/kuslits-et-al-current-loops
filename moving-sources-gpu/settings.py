import shtns
import numpy as np
import torch
from dataclasses import dataclass


@dataclass(frozen=True)
class Geometry:
    r_samples: torch.tensor
    R_oc: float = 3.48e6  # Outer core radius
    r0: float = 0.0  # Normalized inner spherical domain bound.
    r1: float = 1.0  # Normalized outer spherical domain bound.
    nr_r_samples: int = 100  # Number of radial sample (observation) points
    # Geographic resolution (latitude and longitude)
    lat_max: int = 90
    long_max: int = 180


@dataclass(frozen=True)
class PhysicalParams:
    # Permeability of empty space
    mu_0: float = 4 * np.pi * 1e-7  # [H/m]
    # Core conductivity
    sigma: float = 1e6  # older estimate is 5e5 [S/m]
    eta: float = 1 / (mu_0 * sigma)
    # Angular velocity of the rotating flow
    # Omega = 1.2e-10  # FOR AVG OUTER CORE VELOCITIES: v_oc(~10[km/yr])/r_oc(~2500[km])[1/s]
    Omega: float = 5 * 1.2e-10  # FOR MAX POTENTIAL OUTER CORE VELOCITIES: v_oc(~50[km/yr])/r_oc(~2500[km])[1/s]


@dataclass(frozen=True)
class SphericalHarmonicDesc:
    beta_max: int  # Maximum primary source mf two-indexes
    gamma_max: int  # Maximum induced mf two-indexes

    # Spherical harmonic resolution (degree and order)
    lmax: int = 7
    mmax: int = 7

    # First index of the iteration of SH indices
    gamma_first_index: int = 1

    # Alpha indices - for the toroidal flow coefficient
    n: int = 1  # SH degree
    m: int = 0  # SH order


@dataclass(frozen=True)
class ShtnsObjects:
    # Spherical harmonic transformations
    nlat: shtns.sht
    nphi: shtns.sht
    el: shtns.sht
    em: shtns.sht
    l2: shtns.sht


@dataclass(frozen=True)
class GalerkinRepr:
    band_nr: int = 6  # Maximum degree of basis functions
    resolution: int = 100  # Resolution of numerical integrations