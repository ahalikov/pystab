# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

class ChaplyginSledge:
    q = ['phi', 'x', 'y']
    params = ['M', 'R', 'a']
    L = 'M/2*((diff(x(t),t))^2 + (diff(y(t),t))^2 + (R*diff(phi(t),t))^2)'
    Q = {}
    constraints = {'y': 'a*diff(phi(t),t) - diff(x(t),t)*sin(phi(t)) + diff(y(t),t)*cos(phi(t))'}