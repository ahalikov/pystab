# coding=utf-8
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="Дина"
__date__ ="$30.03.2010 22:59:08$"

from pystab.mechanics import *
from sympy import *

def sledge_main():
    sledge = MechanicalFrame('Chaplygin Sledge')

    q, u, a = sledge.add_coordinates('q', 3)
    x, y, phi = q
    dx, dy, dphi = u
    d2x, d2y, d2phi = a
    par = sledge.add_parameters('k, g, alpha')
    k, g, alpha = par
    T = 1/2.* (dx ** 2 + dy ** 2) + 1 / 2. * k ** 2 * dphi ** 2
    sledge.set_vis_viva(T)    
    
    Q = sledge.add_joint_forces({x: g * sin(alpha), y: 0, phi: 0})

    dhc_eqn = [dy - dx * tan(phi)]
    dhc_eqn = solve(dhc_eqn, [dy])
    dhc_eqn = collect(dhc_eqn[dy], [dx, dphi])

    sledge.form_constraints_matrix([dhc_eqn], [dy])
    eqns = sledge.form_voronets_equations()
    print eqns
    
sledge_main()