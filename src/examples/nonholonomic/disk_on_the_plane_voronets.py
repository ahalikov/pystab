# coding=utf-8
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="Дина"
__date__ ="$30.03.2010 23:37:56$"

from pystab.mechanics import *

def disk_main():
    
    disk = MechanicalFrame('Disk on the Plane, p. 256 (301), Markeev')

    q, u, a = disk.add_coordinates('q', 5)
    x, y, phi, psi, theta = q
    dx, dy, dphi, dpsi, dtheta = u
    d2x, d2y, d2phi, d2psi, d2theta = a
    par = disk.add_parameters('m, rho, g')
    m, rho, g = par
    T = 1/2.* m * (dx ** 2 + dy ** 2) + 1 / 8. * m * rho ** 2 * (1 + 4 * cos(theta) ** 2) * dtheta ** 2 + 1 / 8. * m * rho ** 2 * sin(theta) ** 2 * dpsi ** 2 + 1 / 4. * m * rho ** 2 * (dpsi * cos(theta) + dphi) ** 2
    disk.set_vis_viva(T)

    Q = disk.add_joint_forces({x: 0, y: 0, phi: 0, psi: 0, theta: -m * g * rho * cos(theta)})

    dhc_eqn = [dx - rho * (dtheta * sin (psi)  * sin(theta) - (dpsi * cos(theta) + dphi) * cos(psi)), dy + rho * (dtheta * cos(psi) * sin(theta) + (dpsi * cos(theta) + dphi) * sin(psi))]
    dhc_eqn = solve(dhc_eqn, [dx, dy])
    dhc_eqn[dx] = collect(dhc_eqn[dx], [dphi, dpsi, dtheta])
    dhc_eqn[dy] = collect(dhc_eqn[dy], [dphi, dpsi, dtheta])

    disk.form_constraints_matrix([dhc_eqn[dx], dhc_eqn[dy]], [dx, dy])
    eqns = disk.form_voronets_equations()
    print eqns
disk_main()