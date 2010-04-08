# coding=utf-8
__author__="Дина"
__date__ ="$06.04.2010 12:44:25$"
from pystab.deformable import *

def robot_main():
    robot = DeformableFrame('Robot with Two Deformable Wheels', 1, 2)

    q, u, a = robot.define_euler_angles('psi, phi, chi')
    psi, dpdi = [], []
    phi, dphi = [], []
    chi, dchi = [], []
    psi, phi, chi = q['psi'], q['phi'], q['chi']
    dpsi, dphi, dchi = u['psi'], u['phi'], u['chi']

    q, u, a = robot.add_coordinates('q', 2)
    x, y  = q
    dx, dy = u
    d2x, d2y = a

    par = robot.add_parameters('A, l, R1, R2 m1, m2, m, J, C')
    A, l, R1, R2, m1, m2, m, J, C = par
    T =  1/2.* m * ((dx * cos(psi[0]) + dy * sin(psi[0]) + A * dpsi[0]) ** 2 + (-dx * sin(psi[0]) + dy * cos(psi[0])) ** 2) + 1/2. * m1 * ((dx * cos(psi[0]) + dy * sin(psi[0]) - l * dpsi[0]) ** 2 + (-dx * sin(psi[0]) + dy * cos(psi[0])) ** 2) + 1/2. * m2 * ((dx * cos(psi[0]) + dy * sin(psi[0]) + l * dpsi[0]) ** 2 + (-dx * sin(psi[0]) + dy * cos(psi[0])) ** 2) + 1/2. * J * dpsi[0] ** 2 + 1/2. * C * (dphi[0] ** 2 + dphi[1] ** 2)
    robot.set_lagrangian(T)
    
    robot.define_pot_energy()
    MFT = robot.define_main_force_torque()
    def_eqns = robot.define_deformation_eqns([dx * cos(psi[0]) + dy * sin(psi[0]) - cos(chi[0]) * (l * dpsi[0] + R1 * dphi[0]), dx * cos(psi[1]) + dy * sin(psi[1]) + cos(chi[1]) * (l * dpsi[1] - R2 * dphi[1])], [-dx * sin(psi[0]) + dy * cos(psi[0]), -dx * sin(psi[1]) + dy * cos(psi[1])])

    #robot.add_joint_forces({x: 0, psi: 0, y:0, phi0: 0, phi1: 0})
    #eqns = robot.form_lagranges_equations()
    #print 'eqns = ', eqns

robot_main()