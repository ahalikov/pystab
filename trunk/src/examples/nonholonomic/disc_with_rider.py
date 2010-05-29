# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

from pystab.stability import *
from pystab.mechanics import *

disc = MechanicalFrame(name='Disc with rider')

"""
Координаты
theta:
chi:
psi:
phi:
x:
y:
"""
q, u, a = bb.add_coordinates('q', 6)
theta, chi, psi, phi, x, y = q
dtheta, dchi, dpsi, dphi, dx, dy = u
d2theta, d2chi, d2psi, d2phi, d2x, d2y = a

"""
Параметры системы
"""
params = disc.add_parameters('M R m r l mu rho g')
M, R, m, r, l, mu, rho, g = params

"""
Обобщенные силы
"""
Q1 = -M*g*R*sin(theta) - m*g(R+l)*sin(theta) - mu*g*(R+r)*sin(theta) \
    - mu*g*rho*sin(chi - theta)
Q2 = mu*g*rho*sin(chi - theta)

Q = disc.add_joint_forces({theta: Q1, chi: Q2, psi: 0, phi: 0, x: 0, y: 0})

"""
Лагранжиан
"""
L = 1/2*((M*R**2 + m*R*(R+l) + mu*R*(R+r))*(dtheta**2 + dphi**2*cos(theta)**2) \
    + (mu*R*rho)*(dphi*sin(theta) + dpsi)**2) + M/2*(R**2*(dtheta)**2 + 2*R*(dy*cos(phi) \
    - dx*sin(phi))*dtheta*cos(theta) + (dx - R*dphi*sin(theta)*cos(phi))**2 \
    + (dy - R*dphi*sin(theta)*sin(phi))**2) + m/2*((R+l)**2*dtheta**2 + 2*(R+l)*(dy*cos(phi) \
    - dx*sin(phi))*dtheta*cos(theta) + (dx - (R+l)*dphi*sin(theta)*cos(phi))**2 \
    + (dy - (R+l)*dphi*sin(theta)*sin(phi))**2) \
    + mu/2*((R+r)**2*dtheta**2 + rho**2*(dchi - dtheta)**2 \
    + 2*rho*(R+r)*(dchi - dtheta)*dtheta*cos(chi) \
    + 2*(((R+r)*dtheta*cos(theta) + rho*(dchi - dtheta)*cos(chi - theta))*(dy*cos(phi) - dx*sin(phi))) \
    + (dx - dphi*cos(phi)*((R+r)*sin(theta) + rho*sin(chi - theta)))**2 \
    + (dy - dphi*sin(phi)*((R+r)*sin(theta) + rho*sin(chi - theta)))**2)

disc.set_lagrangian(L)

"""
Уравнения связей
"""
dhc_eqns = [dx + R*dpsi*cos(phi), dy + R*dpsi*sin(phi)]


