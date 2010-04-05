# coding=utf-8
from pystab.mechanics import *

def dtwr_main():

    print 'DTWR'
    dtwr = MechanicalFrame('Classic Double Wheeled Robot')

    # before substitute dy = dx*tan(psi)
    """
    q, u, a = dtwr.add_coordinates('q', 5)
    x, y, psi, phi1, phi2  = q
    dx, dy, dpsi, dphi1, dphi2 = u
    d2x, d2y, d2psi, d2phi1, d2phi2 = a
    """
    # after substitute dy = dx*tan(psi)
    q, u, a = dtwr.add_coordinates('q', 4)
    x, psi, phi1, phi2  = q
    dx, dpsi, dphi1, dphi2 = u
    d2x, d2psi, d2phi1, d2phi2 = a
    par = dtwr.add_parameters('A, l, R, m1, m2, m, J, C')
    A, l, R, m1, m2, m, J, C = par
    # Vis-viva before substitute dy = dx*tan(psi)
    #T = 1/2.* m * ((dx * cos(psi) + dy * sin(psi) + A * dpsi) ** 2 + (-dx * sin(psi) + dy * cos(psi)) ** 2) + 1/2. * m1 * ((dx * cos(psi) + dy * sin(psi) - l * dpsi) ** 2 + (-dx * sin(psi) + dy * cos(psi)) ** 2) + 1/2. * m2 * ((dx * cos(psi) + dy * sin(psi) + l * dpsi) ** 2 + (-dx * sin(psi) + dy * cos(psi)) ** 2) + 1/2. * J * dpsi ** 2 + 1/2. * C * (dphi1 ** 2 + dphi2 ** 2)
    # Vis-viva before substitute dy = dx*tan(psi)
    T = 1/2. * m * (dx / cos(psi) + A * dpsi) ** 2 + 1/2. * m1 * (dx / cos(psi) - l * dpsi) ** 2 + 1/2. * m2 * (dx / cos(psi) + l * dpsi) ** 2  + 1/2. * J * dpsi ** 2 + 1/2. * C * (dphi1 ** 2 + dphi2 ** 2)
    dtwr.set_vis_viva(T)

    Q = dtwr.add_joint_forces({x: 0, psi: 0, phi1: 0, phi2: 0})

    # before substitute dy = dx*tan(psi)
    #dhc_eqn = [dx * cos(psi) + dy * sin(psi) - l * dpsi - R * dphi1, dx * cos(psi) + dy * sin(psi) + l * dpsi - R * dphi2]
    # after substitute dy = dx*tan(psi)
    dhc_eqn = [dx * 1 / cos(psi) - l * dpsi - R * dphi1, dx * 1 / cos(psi) + l * dpsi - R * dphi2]
    dhc_eqn = solve(dhc_eqn, [dphi1, dphi2])
    dhc_eqn[dphi1] = collect(dhc_eqn[dphi1], [dx, dpsi])
    dhc_eqn[dphi2] = collect(dhc_eqn[dphi2], [dx, dpsi])

    dtwr.form_constraints_matrix([dhc_eqn[dphi1], dhc_eqn[dphi2]], [dphi1, dphi2])
    eqns = dtwr.form_voronets_equations()
    print eqns

dtwr_main()