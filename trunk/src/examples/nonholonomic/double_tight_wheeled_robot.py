# coding=utf-8
from pystab.mechanics import *

def dtwr_main():

    print 'DTWR'
    dtwr = MechanicalFrame('Classic Double Wheeled Robot')
    
    q, u, a = dtwr.add_coordinates('q', 5)
    x, y, psi, phi1, phi2  = q
    dx, dy, dpsi, dphi1, dphi2 = u
    d2x, d2y, d2psi, d2phi1, d2phi2 = a    
   
    par = dtwr.add_parameters('A, l, R, m1, m2, m, J, C')
    A, l, R, m1, m2, m, J, C = par
    T = 1/2.* m * ((dx * cos(psi) + dy * sin(psi) + A * dpsi) ** 2 + (-dx * sin(psi) + dy * cos(psi)) ** 2) + 1/2. * m1 * ((dx * cos(psi) + dy * sin(psi) - l * dpsi) ** 2 + (-dx * sin(psi) + dy * cos(psi)) ** 2) + 1/2. * m2 * ((dx * cos(psi) + dy * sin(psi) + l * dpsi) ** 2 + (-dx * sin(psi) + dy * cos(psi)) ** 2) + 1/2. * J * dpsi ** 2 + 1/2. * C * (dphi1 ** 2 + dphi2 ** 2)
    dtwr.set_kinetic_energy(T)
    
    Q = dtwr.add_joint_forces({x: 0, psi: 0, y:0, phi1: 0, phi2: 0})

    dhc_eqn = [dx - (l * dpsi + R * dphi1) * cos(psi), dy - (l * dpsi + R * dphi1) * sin(psi), dphi2 - dphi1 - 2 * l * dpsi / R]
    dhc_eqn = solve(dhc_eqn, [dx, dy, dphi2])
    dhc_eqn[dx] = collect(dhc_eqn[dx], [dphi1, dpsi])
    dhc_eqn[dy] = collect(dhc_eqn[dy], [dphi1, dpsi])
    dhc_eqn[dphi2] = collect(dhc_eqn[dphi2], [dphi1, dpsi])

    dtwr.form_constraints_matrix([dhc_eqn[dx], dhc_eqn[dy], dhc_eqn[dphi2]], [dx, dy, dphi2])
    eqns = dtwr.form_voronets_equations()
    #print 'eqns = ', eqns
    
#    eqns['dx'] = dhc_eqn[dx]
#    eqns['dy'] = dhc_eqn[dy]
#    eqns['dphi2'] = dhc_eqn[dphi2]

    q0, u0 = dtwr.define_equilibrium_point(eqns)
    q10, q20, q30, q40, q50 = q0
    u10, u20, u30, u40, u50 = u0
#
#    #test print.to remove
    print q0, u0
#
    manifold = dtwr.form_equilibrium_manifold_equations(eqns)
#
    pertubed = dtwr.form_perturbed_equations(eqns, manifold)

    #print 'Pertubed Eqns: ', pertubed

    fa = dtwr.form_first_approximation_equations(eqns, q0, u0)

    print '1st approximation eqns: ', fa

dtwr_main()