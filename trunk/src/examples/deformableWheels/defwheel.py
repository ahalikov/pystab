# coding=utf-8
# To change this template, choose Tools | Templates
# and open the template in the editor.

from pystab.deformable import *

def wheel_main():
    wheel = DeformableFrame("Wheel", 1, 1)


    q, u, a = wheel.define_euler_angles('psi, phi, chi')
    psi, dpdi = [], []
    phi, dphi = [], []
    chi, dchi = [], []
    psi, phi, chi = q['psi'], q['phi'], q['chi']
    dpsi, dphi, dchi = u['psi'], u['phi'], u['chi']
    print 'psi = ', psi

    q, u, a = wheel.add_coordinates('q', 2)
    x, y  = q
    dx, dy = u

    par = wheel.add_parameters('A, R, m, C, N')
    A, R, m, C, N = par
    T = 1/2.* m * (dx ** 2 + dy ** 2) + 1/ 2. * A * (dpsi[0] ** 2 + dchi[0] **2) + dphi[0] * C * chi[0] * dpsi[0]
    print 'T = ', T
    wheel.set_lagrangian(T)

    wheel.define_pot_energy()
    MFT = wheel.define_main_force_torque()
    def_eqns = wheel.define_deformation_eqns([dx * cos(psi[0]) + dy * sin(psi[0])], [dx * sin(psi[0]) - dy * cos(psi[0])])
    
    wheel.add_joint_forces({x:MFT['Fy'][0] * sin(psi[0]), psi[0]: MFT['Mz'][0], y:-MFT['Fy'][0] * cos(psi[0]), chi[0]: MFT['Mx'][0] - R * MFT['Fy'][0]})
    eqns = wheel.form_lagranges_equations()
    print 'eqns = ', eqns
    
wheel_main()


