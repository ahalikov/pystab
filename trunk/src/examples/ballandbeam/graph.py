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
    [0, 0, 0, 0, 0, 0.1375],
    [0, 0, 0, -45000.0, 0, -5000.0],
    [0, 0, -7.0, 0, 0, 0],
    [-0.759544216220503, 0, 0, 3945.68424010651, 0, 0]
])

B = matrix([0, 0, 0, 1/0.2e-3, 0, 0]).transpose()

u1 = matrix([1.994326034246608, 0.116028092267172, -18.251892618346,
    -0.134385383269116, 2.43091742881793, -0.910486143326368])

u2 = matrix([1.00170584503838, -0.00901093074980716, -20.7834520384088,
    -0.0918151131899239, 2.63748888682019, -0.418876686615918])

def f1(x, t):
    dx = A * matrix(x).transpose()
    return mtx2row(dx)

def f2(x, t):
    dx = (A + B * u2) * matrix(x).transpose()
    return mtx2row(dx)

x0 = [0.03, 0, 0, 0, 0, 0]
slv = scipy_odeint(f2, x0, t0=0, t1=50, last=False, h=1e-2)

#print slv
#print len(slv[:,0])
#print slv[:,1]
t = slv[:, 0]
x1 = slv[:, 1]
x2 = slv[:, 2]
x3 = slv[:, 3]
x4 = slv[:, 4]
dx1 = slv[:, 5]
dx2 = slv[:, 6]

plt.figure(1)
plt.subplot(221)
plt.plot(t, x1, label=r'Ball position $\rho$')
plt.xlabel('t')
plt.ylabel(r'$\rho$')
plt.grid(True)
plt.legend(loc=2)

plt.subplot(222)
plt.plot(t, x2, label=r'Gear angle $\theta$', color='green')
plt.plot(t, x3, label=r'Beam angle $\alpha$', color='blue')
plt.xlabel('t')
plt.grid(True)
plt.legend(loc=2)

plt.subplot(223)
plt.plot(t, x4, label=r'Current $\gamma$', color='red')
plt.xlabel('t')
plt.ylabel(r'$\gamma$')
plt.grid(True)
plt.legend(loc=3)

plt.subplot(224)
plt.plot(t, dx1, label=r'Ball velocity $\dot_\rho$')
plt.plot(t, dx2, label=r'Gear velocity $\dot_\theta$', color='green')
plt.xlabel('t')
plt.grid(True)
plt.legend(loc=4)

#t, x3, 'blue', t, x4, 'orange', t, dx1, 'gray', t, dx2, 'black'

#plt.plot(t, x1, 'red', t, x2, 'green', t, x3, 'blue', t, x4, 'orange', t, dx1, 'gray', t, dx2, 'black')
plt.grid(True)
plt.show()

#End