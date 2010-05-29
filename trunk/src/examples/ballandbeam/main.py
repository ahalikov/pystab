# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="Artur"
__date__ ="$03.02.2010 12:55:33$"

"""
from numpy import matrix
from pystab.stability import *
from pystab.integration.ode import *
import matplotlib.pyplot as plt
"""

def main():
    #import examples.ballandbeam.graph

    #import examples.ballandbeam.ballandbeam
    import examples.nonholonomic.chaplygin_sledge_voronets

    """
    A = matrix([[0, 1, 0], [1, 1, 1], [0, 0, 0]])
    B = matrix([0, 1, 1]).transpose()

    lq = LQRegulator(A, B)
    u = lq.find_control()
    print u

    def f(x, t):
        dx = (A  + B * u) * matrix(x).transpose()
        return mtx2row(dx)

    x0 = [1e-3, 1e-3, 1e-3]
    slv = scipy_odeint(f, x0, t0=0, t1=20, last=False, h=1e-2)

    t = slv[:, 0]
    x1 = slv[:, 1]
    x2 = slv[:, 2]
    x3 = slv[:, 3]

    plt.figure(1)
    plt.plot(t, x1, 'red', t, x2, 'green', t, x3, 'blue')
    plt.grid(True)
    plt.show()
    """
    
if __name__ == "__main__":
    main()
