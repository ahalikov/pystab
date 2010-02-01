# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

from sympy.matrices.matrices import eye
from sympy import Matrix, matrix2numpy, zeros
from numpy import sum, where, linspace, array
from numpy.linalg import svd
from scipy import integrate

def ctrb(A, B):
    """
    Calculates the controllability matrix for pair A, B.
    """
    n = A.shape[0]
    m = B.shape[1]
    assert n == A.shape[1] and n == B.shape[0]
    C = zeros([n, n])
    C[:, 0] = B
    i = 1
    while i < n:
        C[:, i] = A * C[:, i-1]
        i += 1
    return C

def matrix_rank(A, tol=1e-8):
    s = svd(matrix2numpy(A), compute_uv=0)
    return sum(where(s > tol, 1, 0))

def row2mtx(row, n, debug=False):
    mtx = [row[i*n : (i+1)*n] for i in range(n)]
    if debug: print mtx
    return Matrix(mtx)

def mtx2row(mtx, numpyarr=False):
    row = [i for i in mtx]
    if numpyarr:
        row = array(row)
    return row

def is_controllable(A, B):
    rank = matrix_rank(ctrb(A, B))
    return A.shape[0] == rank

class LQRegulator:

    def __init__(self, A, B, time=5):

        assert A.is_square == True
        assert A.shape[0] == B.shape[0]

        self.sys_dim = A.shape[0]
        self.control_dim = B.shape[1]
        self.A = A
        self.B = B
        self.alpha = eye(self.sys_dim)
        self.beta = eye(self.control_dim)
        self.control = []

    def find_control(self, time=5):
        """
        Calculates the optimal control.
        """
        n = self.B.shape[0]
        # Checking controllability
        assert is_controllable(self.A, self.B) == True
        t = linspace(0, time, 10000)
        # Initial values
        c0 = array([0 for i in range(n*n)])
        # Integration of diff. equations
        y = integrate.odeint(self.deriv, c0, t)
        # Take only last values
        y = row2mtx(y[y.shape[0]-1].tolist(), n)
        self.C = y
        self.control = -self.beta.inv() * self.B.transpose() * y
        return self.control

    def deriv(self, C, t):
        """
        Right part of the Lyapunov-Bellman-Riccati equation
        """
        m = self.B.shape[0]
        BBiBt = self.B * self.beta.inv() * self.B.transpose()
        C_mtx = row2mtx(Matrix(C), m)
        res = -C_mtx*BBiBt*C_mtx + self.A.transpose()*C_mtx + C_mtx*self.A + self.alpha
        return mtx2row(res)

    def get_control(self):
        return self.control

#End