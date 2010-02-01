# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

class DiscWithRider:
    q = ['theta', 'chi', 'psi', 'phi', 'x', 'y']

    params = ['M', 'R', 'm', 'r', 'l', 'mu', 'rho', 'g']

    Q = {
    'theta': '-M*g*R*sin(theta(t)) - m*g(R+l)*sin(theta(t))\
    - mu*g*(R+r)*sin(theta(t)) - mu*g*rho*sin(chi(t) - theta(t))',
    'chi': 'mu*g*rho*sin(chi(t)-theta(t))',
    'psi': 0,
    'phi': 0,
    'x': 0, 'y': 0}

    L = '1/2*((M*R^2+m*R*(R+l)+mu*R*(R+r))*(diff(theta(t),t)^2\
    + diff(phi(t),t)^2*cos(theta(t))^2) + (mu*R*rho)*(diff(phi(t),t)*sin(theta(t))\
    + diff(psi(t),t))^2) + M/2*(R^2*(diff(theta(t),t))^2 + 2*R*(diff(y(t),t)*cos(phi(t))\
    - diff(x(t),t)*sin(phi(t)))*diff(theta(t),t)*cos(theta(t))\
    + (diff(x(t),t) - R*diff(phi(t),t)*sin(theta(t))*cos(phi(t)))^2\
    + (diff(y(t),t) - R*diff(phi(t),t)*sin(theta(t))*sin(phi(t)))^2)\
    + m/2*((R+l)^2*diff(theta(t),t)^2 + 2*(R+l)*(diff(y(t),t)*cos(phi(t))\
    - diff(x(t),t)*sin(phi(t)))*diff(theta(t),t)*cos(theta(t))\
    + (diff(x(t),t) - (R+l)*diff(phi(t),t)*sin(theta(t))*cos(phi(t)))^2\
    + (diff(y(t),t) - (R+l)*diff(phi(t),t)*sin(theta(t))*sin(phi(t)))^2)\
    + mu/2*((R+r)^2*diff(theta(t),t)^2 + rho^2*(diff(chi(t),t) - diff(theta(t),t))^2\
    + 2*rho*(R+r)*(diff(chi(t),t) - diff(theta(t),t))*diff(theta(t),t)*cos(chi(t))\
    + 2*(((R+r)*diff(theta(t),t)*cos(theta(t)) + rho*(diff(chi(t),t)\
    - diff(theta(t),t))*cos(chi(t) - theta(t)))*(diff(y(t),t)*cos(phi(t))\
    - diff(x(t),t)*sin(phi(t)))) + (diff(x(t),t)\
    - diff(phi(t),t)*cos(phi(t))*((R+r)*sin(theta(t)) + rho*sin(chi(t) - theta(t))))^2\
    + (diff(y(t),t) - diff(phi(t),t)*sin(phi(t))*((R+r)*sin(theta(t)) + rho*sin(chi(t) - theta(t))))^2)'

    constraints = {
        'x': 'diff(x(t),t) + R*diff(psi(t),t)*cos(phi(t))',
        'y': 'diff(y(t),t) + R*diff(psi(t),t)*sin(phi(t))'
    }

