# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

from numpy import array, matrix, zeros, eye, sum, where
from numpy.linalg import svd, inv
from pystab.integration.ode import *
import sympy

def ctrb(A, B):
    """
    Calculates the matrix of controllability for pair A, B.
    """
    n = A.shape[0]
    assert n == A.shape[1] and n == B.shape[0]
    C = matrix(zeros([n, n]))
    C[:, 0] = B
    i = 1
    while i < n:
        C[:, i] = A * C[:, i-1]
        i += 1
    return C

def sym_ctrb(A, B):
    """
    Calculates the matrix of controllability for pair A, B in symbolic form.
    """
    n = A.shape[0]
    assert n == A.shape[1] and n == B.shape[0]
    C = sympy.zeros([n, n])
    C[:, 0] = B
    i = 1
    while i < n:
        C[:, i] = A * C[:, i-1]
        i += 1
    return C

def matrix_rank(A, tol=1e-8):
    """
    Calculates the rank of a matrix.
    """
    s = svd(A, compute_uv=0)
    return sum(where(s > tol, 1, 0))

def row2mtx(row, n):
    """
    Row-to-matrix conversion.
    [a, b, c, d] -> [[a, b], [c, d]]
    """
    mtx = matrix(zeros([n, n]))
    for i in range(n):
        for j in range(n):
            if isinstance(row, matrix):
                mtx[i, j] = row[0, i*n + j]
            else:
                mtx[i, j] = row[i*n + j]
    return mtx

def mtx2row(mtx):
    """
    Matrix-to-row conversion.
    [[a, b], [c, d]] -> [a, b, c, d]
    """
    row = []
    for i in range(mtx.shape[0]):
        for j in range(mtx.shape[1]):
            row.append(mtx[i, j])
    return row

def is_controllable(A, B):
    """
    Checks the controllability of pair A, B.
    """
    rank = matrix_rank(ctrb(A, B))
    return A.shape[0] == rank

class LQRegulator:
    """
    Linear - quadratic regulator.
    """
    def __init__(self, A, B):

        assert A.shape[0] == A.shape[1]
        assert A.shape[0] == B.shape[0]

        self.sys_dim = A.shape[0]
        self.control_dim = B.shape[1]
        self.A = A
        self.B = B
        self.alpha = eye(self.sys_dim)
        self.beta = eye(self.control_dim)
        self.control = []

    def find_control(self, time=10, h=1e-3):
        """
        Calculates the optimal control.
        """
        # Checking controllability
        #if not is_controllable(self.A, self.B):
        #    raise ArgumentError, "Pair A,B is not controllable."
        # Initial values
        C0 = array([0 for i in range(self.sys_dim**2)])
        # Integration of diff. equations
        C = scipy_odeint(self.deriv, C0, time, h)
        self.C = row2mtx(C, self.sys_dim)
        self.control = -inv(self.beta) * self.B.transpose() * self.C
        return self.control

    def deriv(self, C, t):
        """
        Right part of the Lyapunov-Bellman-Riccati equation.
        """
        BBiBt = self.B * inv(self.beta) * self.B.transpose()
        C_mtx = row2mtx(matrix(C), self.sys_dim)
        res = -C_mtx*BBiBt*C_mtx + self.A.transpose()*C_mtx + C_mtx*self.A + self.alpha
        return mtx2row(res)

    def get_control(self):
        return self.control

"""
For quick tests
"""
def main():
    A = matrix([[0, 1, 0], [1, 0, 1], [0, 0, 0]])
    B = matrix([[0, 0], [0, 1], [1, 0]])
    C  = ctrb(A, B)
    print C
    #lq = LQRegulator(A, B)
    #print lq.find_control()

if __name__ == '__main__':
    main()

#End