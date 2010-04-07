# coding=utf-8
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="Дина"
__date__ ="$05.04.2010 11:03:34$"

from mechanics import *
from ctypes import ArgumentError
from sympy import (Symbol, Function, symbols, sin, cos, tan, cot,
    solve, zeros, Derivative as D, diff, Eq, collect, Matrix, pprint,
    simplify, radsimp, fraction, together, expand)
import const


t = Symbol('t')

# constants
const.Simple_Keldysh = 1
const.Levin_Fufaev = 2
const.Big_Keldysh = 3
const.Simple_Keldysh_Consts = 'alpha, beta, gamma, N, cy, ct, nu2, rho2'
const.Levin_Fufaev_Consts = ''
const.Big_Keldysh_Consts = ''
const.Force_Torque_List = 'Fx, Fy, Fz, Mx, My, Mz, '
const.Euler_Angles = 'psi, phi, chi'


class DeformableFrame (MechanicalFrame):

    def __init__(self, name='', theory = const.Simple_Keldysh, wheelCnt = 1):
        '''
        Theory parameter takes integer values
        1 - takes into account only lateral deformation
            (according to article written by Neymark, Fufaev in 1971)
        2 - takes into account both lateral and longitudinal deformations
            (according to chapter 1 of book written by Levin, Fufaev in 1989)
        3 - takes into account both lateral and longitudinal deformations
            (according to chapter 4.2 of book written by Levin, Fufaev in 1989)

        WheelCnt parameter takes integer for deformable wheels count
        '''
        MechanicalFrame.__init__(self, name)
        self.theory = theory
        self.wheels_count = wheelCnt

        self.main_force_torque = {}

        self.euler_q = {}
        self.euler_u = {}
        self.euler_a = {}
        self.def_const = {}
        unchanged_params = ''
        self.deformation_params = ''

        self.deformation_q = {}
        self.deformation_u = {}
        self.deformation_eqns = {}

        # set deformation parameters
        self.__set_deformation_parameters()

        # add unchanged parameters for theory 1
        if self.theory == const.Simple_Keldysh:
            unchanged_params = const.Simple_Keldysh_Consts
            self.deformation_params = 'eta0, eta1, '
        elif self.theory == const.SimpleLevin_Fufaev:
            unchanged_params = const.Levin_Fufaev_Consts
            self.deformation_params = 'eta0, eta1, xi0'
        elif self.theory == const.Big_Keldysh:
            unchanged_params = const.Big_Keldysh_Consts
            self.deformation_params = 'eta0, eta1, xi0, xi1'
        params = list(symbols(unchanged_params))
        print 'self.deformation_params = ', self.deformation_params
        print 'self.theory = ', self.theory
        unchanged_params = ''
        for p in params:
            self.def_const[str(p)] = []
            for i in range(self.wheels_count):
                self.def_const[str(p)].append(Symbol(str(p) + str(i)))
                unchanged_params += str(p) + str(i)
                unchanged_params += ', '
        print 'def_const = ', self.def_const
        self.add_parameters(unchanged_params)

        # set deformation parameters
        self.__set_deformation_parameters()

    def define_main_force_torque(self, string = const.Force_Torque_List):
        '''
        Get string with 6 elements. Defines names for main force and torque
        '''
        params = list(symbols(string))
        for p in params:
            self.main_force_torque[str(p)] = []
            for i in range (self.wheels_count):
                if str(p) == 'Fx':
                    self.main_force_torque[str(p)].append(0)
                if str(p) == 'Fy':
                    self.main_force_torque[str(p)].append(pdiff(self.pot_energy, self.deformation_q['eta0'][i]))
                if str(p) == 'Fz':
                    self.main_force_torque[str(p)].append(0)
                if str(p) == 'Mx':
                    self.main_force_torque[str(p)].append(- pdiff(self.pot_energy, self.euler_q['chi'][i]))
                if str(p) == 'My':
                    self.main_force_torque[str(p)].append(0)
                if str(p) == 'Mz':
                    self.main_force_torque[str(p)].append(pdiff(self.pot_energy, self.deformation_q['eta1'][i]))
        print 'main_force_torque = ', self.main_force_torque

        return self.main_force_torque

    def define_euler_angles(self, string = const.Euler_Angles):
        '''
        Defines Euler Angles for each deformable wheel.
        '''
        params = list(symbols(string))
        for p in params:
            self.euler_q[str(p)] = []
            self.euler_u[str(p)] = []
            self.euler_a[str(p)] = []
            q, u, a = self.add_coordinates('q', self.wheels_count)
            if self.wheels_count > 1:
                for i in range (self.wheels_count):
                    self.euler_q[str(p)].append(q[i])
                    self.euler_u[str(p)].append(u[i])
                    self.euler_a[str(p)].append(a[i])
            else:
                self.euler_q[str(p)].append(q)
                self.euler_u[str(p)].append(u)
                self.euler_a[str(p)].append(a)

        print 'self.euler_q = ', self.euler_q
        print 'self.euler_u = ', self.euler_u
        print 'self.euler_a = ', self.euler_a

        return self.euler_q, self.euler_u, self.euler_a

    def define_pot_energy(self):
        # define pot_energy for each theory
        self.pot_energy = 0
        for i in range(self.wheels_count):
            self.pot_energy += 1./2. * (self.def_const['cy'][i] * self.deformation_q['eta0'][i] ** 2 -\
                2 * self.def_const['nu2'][i] * self.def_const['N'][i] * self.deformation_q['eta0'][i] * \
                self.euler_q['chi'][i] + self.def_const['rho2'][i] * self.def_const['N'][i] * self.euler_q['chi'][i] ** 2\
                + self.def_const['ct'][i] * self.deformation_q['eta1'][i] ** 2)
        print 'P = ', self.pot_energy
        return self.pot_energy


    def __set_deformation_parameters(self):
        '''
        Define deformable parameters that depend on time
        '''
        print 'self.deformation_params = ', self.deformation_params
        def_params = list(symbols(self.deformation_params))
        for p in def_params:
            self.deformation_q[str(p)] = []
            self.deformation_u[str(p)] = []
            q, u, a = self.add_coordinates('q', self.wheels_count)
            if self.wheels_count > 1:
                for i in range(self.wheels_count):
                    self.deformation_q[str(p)].append(q[i])
                    self.deformation_u[str(p)].append(u[i])
            else:
                 self.deformation_q[str(p)].append(q)
                 self.deformation_u[str(p)].append(u)

        print 'self.deformation_q = ', self.deformation_q
        print 'self.deformation_u = ', self.deformation_u

    def define_deformation_eqns(self, vx, vy):
        '''
        Define deformation eqns.
        vx - list of deformable wheels longitudinal velocity
        vy - list of deformable wheels lateral velocity
        '''
        def_params = list(symbols(self.deformation_params))
        for p in def_params:
            self.deformation_eqns[str(p)] = []
        if self.theory == const.Simple_Keldysh:
            for i in range(self.wheels_count):
                self.deformation_eqns['eta0'].append(vy[i] + self.deformation_u['eta0'][i] + \
                        vx[i] * (self.euler_q['psi'][i] + self.deformation_q['eta1'][i]))
                self.deformation_eqns['eta1'].append(self.euler_u['psi'][i] + \
                    self.deformation_u['eta1'][i] + \
                        vx[i] * (self.def_const['alpha'][i] * self.deformation_q['eta0'][i]\
                        + self.def_const['beta'][i] * self.deformation_q['eta1'][i]\
                        + self.def_const['gamma'][i] * self.euler_q['chi'][i]))
        elif self.theory == const.Levin_Fufaev:
            self.deformation_eqns = {}
        elif self.theory == const.Big_Keldysh:
            self.deformation_eqns = {}
        print 'deformation_eqns = ', self.deformation_eqns
        return self.deformation_eqns

def robot_main():
    robot = DeformableFrame('Budanov', 1, 1)
    q, u, a = robot.define_euler_angles()
    print 'q = ', q
    print 'u = ', u
    print 'a = ', a
    robot.define_pot_energy()
    mainFT = robot.add_main_force_torque()
    print 'main force and torque = ', mainFT
