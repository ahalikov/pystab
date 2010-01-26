# coding=windows-1251
__author__="akhalikoff at gmail dot com"
__date__ ="$23.01.2010 9:59:42$"

from ctypes import ArgumentError
from numpy.core.defmatrix import matrix
from sympy.matrices.matrices import eye

def is_controllable(A, B):
    """
    Checks the controllability of pair matrixes A, B
    """
    pass

def ctrb(A, B):
    """
    Calculates the controllability matrix for pair matrixes A, B
    """
    n = A.shape[0]
    m = B.shape[1]
    if n != A.shape[1] or m != B.shape[0]:
        raise ArgumentError, 'Matrix A must be square with as many rows as B.'

    w = zeros([n, n * m]);
    #Matrix tmp = B;
    #for (int k = 0; k < n; k++) {
    #    for (int i = 0; i < n; i++) {
    #        for (int j = m*k; j < m*(k+1); j++) {
    #            W.set(i, j, tmp.get(i, j - m*k));
    #        }
    #    }
    #    tmp = A.times(tmp);
    #}
    return w

class LQRegulator:

    def __init__(self, A, B):
        self.A = A
        self.B = B
        self.m = A.shape[0]
        self.n = B.shape[1]
        self.alpha = eye(m)
        self.beta = eye(n)
        
    