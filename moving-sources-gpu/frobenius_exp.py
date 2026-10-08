import numpy as np
from bessel_settings import Omega, eta, elsasser_gamma
from bessel_settings import n_gamma, p, rho, r

lambda_gamma = 0.1

A_gamma = Omega*elsasser_gamma/(4*np.pi*eta)
B_gamma = lambda_gamma/eta

C1 = 0
C2 = - lambda_gamma/(2*eta*(2*n_gamma+3))
C3 = - Omega*elsasser_gamma/(12*eta*np.pi*(2*n_gamma+4))
C4 = lambda_gamma**2/(8*eta**2*(2*n_gamma+3)*(2*n_gamma+5))
C5 = - (A_gamma*C2+B_gamma*C3)/(5*(2*n_gamma+6))


def frobenius(r):
    frob_n = r**n_gamma*(C2*r**2+C3*r**3+C4*r**4)
    return frob_n

r_nondim = np.linspace(0,1,100)

frob_vals = frobenius(r_nondim)

#print(frob_vals)


def frob_lambda(lambda_gamma):
    r = 0.9
    A_gamma = Omega*elsasser_gamma/(4*np.pi*eta)
    B_gamma = lambda_gamma/eta
    #C1 = 0
    C2 = - lambda_gamma/(2*eta*(2*n_gamma+3))
    C3 = - Omega*elsasser_gamma/(12*eta*np.pi*(2*n_gamma+4))
    C4 = lambda_gamma**2/(8*eta**2*(2*n_gamma+3)*(2*n_gamma+5))
    C5 = - (A_gamma*C2+B_gamma*C3)/(5*(2*n_gamma+6))
    frob_n = r ** n_gamma * (C2 * r ** 2 + C3 * r ** 3 + C4 * r ** 4 + C5 * r**5)
    return frob_n

def frob_prime_lambda(lambda_gamma):
    r = 0.9
    A_gamma = Omega*elsasser_gamma/(4*np.pi*eta)
    B_gamma = lambda_gamma/eta
    C0 = 1
    #C1 = 0
    C2 = - lambda_gamma/(2*eta*(2*n_gamma+3))
    C3 = - Omega*elsasser_gamma/(12*eta*np.pi*(2*n_gamma+4))
    C4 = lambda_gamma**2/(8*eta**2*(2*n_gamma+3)*(2*n_gamma+5))
    C5 = - (A_gamma*C2+B_gamma*C3)/(5*(2*n_gamma+6))
    frob_prime = r ** n_gamma * (C0 +
                                 (2+n_gamma)*C2 * r ** (2-1) +
                                 (3+n_gamma)*C3 * r ** (3-1) +
                                 (4+n_gamma)*C4 * r ** (4-1) +
                                 (5+n_gamma)*C5 * r ** (5-1))
    return frob_prime


def robin_boundary_condition(lambda_gamma):
    return n_gamma*frob_lambda(lambda_gamma) + frob_prime_lambda(lambda_gamma)


lambda_vals = np.linspace(1,100,100)
bc_plot = robin_boundary_condition(lambda_vals)

# import matplotlib.pyplot as plt
# fig1 = plt.figure(1)
# plt.plot(lambda_vals, bc_plot)
# plt.savefig('frobenius_bc_lambda.png')


roots = (np.abs(bc_plot) < max(abs(bc_plot))/5e1)

print(lambda_vals[roots])

def frob_lambda_r(r, lambda_gamma):
    A_gamma = Omega*elsasser_gamma/(4*np.pi*eta)
    B_gamma = lambda_gamma/eta
    C0 = 1
    #C1 = 0
    C2 = - lambda_gamma/(2*eta*(2*n_gamma+3))
    C3 = - Omega*elsasser_gamma/(12*eta*np.pi*(2*n_gamma+4))
    C4 = lambda_gamma**2/(8*eta**2*(2*n_gamma+3)*(2*n_gamma+5))
    C5 = - (A_gamma*C2+B_gamma*C3)/(5*(2*n_gamma+6))
    frob_n = r ** n_gamma * (C0 + C2 * r ** 2 + C3 * r ** 3 + C4 * r ** 4 + C5 * r**5)
    return frob_n

frob_lambda_r_plot = frob_lambda_r(r_nondim, np.array([0]))#lambda_vals[roots])

import matplotlib.pyplot as plt
fig2 = plt.figure(2)
plt.plot(r_nondim, frob_lambda_r_plot)
plt.savefig('frobenius_function_r.png')

print('Mo. vege.')



