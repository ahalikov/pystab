# coding=windows-1251

__author__="Artur"
__date__ ="$21.12.2009 21:36:35$"

from MechanicalFrame import *

def sp_main():

    sp = MechanicalFrame()
    q, u, a = sp.add_coordinates('q')

    phi, = q
    dphi, = u
    d2phi, = a

    p = sp.add_parameters('m,l,r,Ic,g')
    m,l,r,Ic,g = p

    T = 1/4.0*m*(r**2 + 2*(l + r)**2)*dphi**2
    P = -m*g*(l + r)*cos(phi)
    L = T - P

    sp.set_lagrangian(L)

    eqns = sp.form_lagranges_equations()
    pprint(eqns)
    