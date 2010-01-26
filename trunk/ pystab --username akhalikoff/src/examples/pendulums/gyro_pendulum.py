# coding=windows-1251

__author__="akhalikoff at gmail dot com"
__date__ ="$19.01.2010 10:48:37$"

"""
Гироскопический маятник
Введение в аналитическую механику, Н.В.Бутенин, Н.А.Фуфаев, 1991, стр.137.
"""

from MechanicalFrame import *

def gp_main():
    
    gp = MechanicalFrame('Gyro pendulum')
    
    q, u, a = gp.add_coordinates('q', 3)
    theta, psi, phi = q
    dtheta, dpsi, dphi = u
    d2theta, d2psi, d2phi = a

    A, Iy, B, m1, l1, m2, l2, Ix, g = gp.add_parameters('A Iy B m1 l1 m2 l2 Ix g')
    #A = m1*l1^2 + m2*l2^2 + Ix
    #B = (m1*l1 - m2*l2)*g

    T = 1/2.0*(A*(dtheta**2 + (dpsi**2)*sin(theta)**2) + Iy*(dpsi*cos(theta) + dphi)**2)
    P = B*cos(theta)
    L = T - P
    gp.set_lagrangian(L)
    #pprint(L)

    eqns = gp.form_lagranges_equations()
    #pprint(eqns[d2phi])

    q0, u0 = gp.define_equilibrium_point(eqns)
    q10, q20, q30 = q0
    u10, u20, u30 = u0
    #print q0, u0

    manifold = gp.form_equilibrium_manifold(eqns)
    pprint(manifold)

