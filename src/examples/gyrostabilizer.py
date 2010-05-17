from sympy.galgebra.latex_ex import LaTeX
import sympy
# coding=utf-8

__author__="Artur"
__date__ ="$16.05.2010 21:14:57$"

from pystab.mechanics import normalize
from numpy.core.setup import sym2def
from pystab.mechanics import *

gs = MechanicalFrame(name='GyroStabilizer')

# Coordinates
q, u, a = gs.add_coordinates('q', 2)
alpha, beta = q
dalpha, dbeta = u
d2alpha, d2beta = a

"""
System parameters

Inertia:
A1 = Ar + Aj + Aa
A2 = Aa + Cj
A3 = Ar + Aj - Cj
B1 = Br + Bj

Moments:
M1, M2, M3

Velocity of the rotor is constant:
domega = Omega
"""
p = gs.add_parameters('Cr A1 A2 A3 B1, M1, M2, M3 Omega s u')
Cr, A1, A2, A3, B1, M1, M2, M3, Omega, s, u = p

# Motion equations (without control)
meqns = {
    # alpha
    d2alpha: d2alpha*(A1*cos(beta)**2 + A2*sin(beta)**2) -\
        A3*dalpha*dbeta*sin(2*beta) + Cr*(Omega + dalpha*sin(beta))*dbeta*cos(beta),
    # beta
    d2beta: B1*d2beta + A3*dalpha**2*sin(beta)*cos(beta) -\
        Cr*(Omega + dalpha*sin(beta))*dalpha*cos(beta)
    }
#print meqns
#pprint(meqns)

# Normalization
meqns = normalize(meqns)
#pprint(meqns)

q0, u0 = gs.define_equilibrium_point(meqns)
#print q0
#print u0
q0[alpha] = q0[beta] = 0
u0[dalpha] = u0[dbeta] = 0

manifold = gs.form_equilibrium_manifold_equations(meqns, {u:0})
#pprint(manifold)

# Equations of perturbed motion
peqns = gs.form_perturbed_equations(meqns, manifold)
#pprint(peqns)

# Approximated equations of perturbed motion
fa_eqns = gs.form_first_approximation_equations(peqns, q0, simplified=False)
pprint(fa_eqns)