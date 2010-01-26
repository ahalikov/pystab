# coding=utf-8

__author__="Artur"
__date__ ="$09.12.2009 13:27:17$"

class ChaplyginSledge:
    q = ['phi', 'x', 'y']
    params = ['M', 'R', 'a']
    L = 'M/2*((diff(x(t),t))^2 + (diff(y(t),t))^2 + (R*diff(phi(t),t))^2)'
    Q = {}
    constraints = {'y': 'a*diff(phi(t),t) - diff(x(t),t)*sin(phi(t)) + diff(y(t),t)*cos(phi(t))'}