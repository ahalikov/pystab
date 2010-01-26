# coding=windows-1251
__author__="akhalikoff at gmail dot com"
__date__ ="$23.01.2010 9:59:42$"

from numpy.core.defmatrix import matrix
from sympy.matrices.matrices import eye



class LQRegulator:
    def __init__(self, A, B):
        self.A = A
        self.B = B
        self.alpha = eye(3)
        self.beta = eye(2)
