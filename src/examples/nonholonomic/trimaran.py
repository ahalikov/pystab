from pystab.mechanics import *
from sympy import *

def trim_main():

    trim = MechanicalFrame('Trimaran Model')

    q, u, a = trim.add_coordinates('q', 4)
    x, y, phi, psi = q
    dx, dy, dphi, dpsi = u
    d2x, d2y, d2phi, d2psi = a
    par = trim.add_parameters('R, m1, m2, rho, M')
    R, m1, m2, rho, M = par
    T = 1/2.* (m1 + 2 * m2) * (dx ** 2 + dy ** 2) + 1/2. * (2 * m2 * rho ** 2 + m1 * R ** 2 / 2) * dpsi ** 2 + m1 * R ** 2 * dphi ** 2 / 4
    trim.set_vis_viva(T)
    
    Q = trim.add_joint_forces({x: 0, y: 0, phi: M, psi: 0})
    
    dhc_eqn = [dx - R*dphi*cos(psi), dy - R * dphi * sin(psi)]
    dhc_eqn = solve(dhc_eqn, [dx, dy])
    dhc_eqn[dx] = collect(dhc_eqn[dx], [dphi, dpsi])
    dhc_eqn[dy] = collect(dhc_eqn[dy], [dphi, dpsi])   

    trim.form_constraints_matrix([dhc_eqn[dx], dhc_eqn[dy]], [dx, dy])
    eqns = trim.form_voronets_equations()
    print eqns
trim_main()