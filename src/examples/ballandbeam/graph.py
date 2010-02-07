# coding=utf-8

__author__="artur"
__date__ ="$01.02.2010 19:05:06$"

from numpy import matrix
import matplotlib.pyplot as plt
from pylab import *

from pystab.integration.ode import *
from pystab.stability import *

A = matrix([
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0.137500000000000],
    [0, 0, 0, -45000.0000000000, 0, -500.000000000000],
    [0, 0, -7.00000000000000, 0, 0, 0],
    [-0.759544216220503, 0, 0, 3945.68424010651, 0, 0]
])

B = matrix([0, 0, 0, 1/0.2e-3, 0, 0]).transpose()

u = matrix([1.994326034246608, 0.116028092267172, -18.251892618346,
    -0.134385383269116, 2.43091742881793, -0.910486143326368])

def f1(x, t):
    dx = A * matrix(x).transpose()
    return mtx2row(dx)

def f2(x, t):
    dx = A * matrix(x).transpose() + B * u
    return mtx2row(dx)

x0 = [1e-3, 0.004, 2e-3, 1e-4, 1e-4, 1e-4]
slv = scipy_odeint(f2, x0, t0=0, t1=5, last=False, h=1e-1)

print slv

#plt.figure(1)
#plt.plot(slv[:,0], slv[:,5])
#plt.grid(True)
#plt.show()

#End