# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

from pystab.mechanics import *

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
    