# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

"""
Ball on a Beam system.
Система "Шарик на желобе".
"""

from pystab.stability import matrix_rank
import numpy
from numpy import linalg as la
import matplotlib.pyplot as plt

from pystab.stability import *
from pystab.mechanics import *

bb = MechanicalFrame(name='Ball on a Beam')

"""
Обобщенные координаты
rho - положение шарика на желобе, отсчитываемое от левого неподвижного
конца желоба;
theta - угловая скорость колеса двигателя;
alpha - угол на который поднимается желоб.
"""
q, u, a = bb.add_coordinates('q', 3)
rho, theta, alpha = q
drho, dtheta, dalpha = u
d2rho, d2theta, d2alpha = a

"""
Параметры системы:

m, r, J - масса, радиус и момент инерции шарика
M, l - масса и длина желоба
l1 - длина рычага, соединяющего колесо двигателя с желобом
R, Jm - радиус и момент инерции колеса
g - гравитационная постоянная
k = (5 + 2*(R/R')^2)/5, где R' - малый радиус качения шарика.
Kg - передаточное число
rho0 - исследуемое положение шарика
tau - управляющий момент на двигателе (tau = Kg * K2 * i)
где K2 - это постоянная момента двигателя, i - ток.
"""
p1 = bb.add_parameters('m r J M l l1 R Jm g k Kg rho0 tau')
m, r, J, M, l, l1, R, Jm, g, k, Kg, rho0, tau = p1

"""
Q - обобщенные силы,
T - кинетическая энергия,
P - потенциальная энергия,
L - функция Лагранжа.
"""
Q = bb.add_joint_forces({r: 0, theta: tau, alpha: 0})
T = 1.0/2*(k*m*drho**2 + (J + m*(rho**2))*dalpha**2 + Jm*dtheta**2)
#printm(T)
P = g*sin(alpha)*(m*rho + 1.0/2*M*l)
#printm(P)
L =  T - P
#printm(L)
bb.set_lagrangian(L)

"""
Полное и упрощенное (alpha~0, theta~0) уравнения связи.
"""
dhc_eqn_full = [l**2*dalpha*sin(alpha) + R**2*dtheta*sin(theta)+ \
    l*R*(dalpha - dtheta)*sin(alpha - theta) - l*R*(dtheta*sin(theta)+ \
    dalpha*sin(alpha)) + l*l1*dalpha*cos(alpha) - l1*R*dtheta*cos(theta)]

dhc_eqn_simple = [l*dalpha - R*dtheta]
dhc_eqn = dhc_eqn_full
#printm(dhc_eqn)

"""
Уравнение связи, разрешенное относительно зависимой скорости.
"""
dhc_eqn = solve(dhc_eqn, [dalpha])
dhc_eqn = collect(dhc_eqn[dalpha], dtheta)
bb.form_constraints_matrix([dhc_eqn], [dalpha])

"""
Получаем уравнения Шульгина для r, theta и alpha.
Флаг first_order означает, что уравнения будут сразу приведены к первому порядку,
и как следствие, увеличится размерность системы.
"""
eqns = bb.form_shulgins_equations(normalized=True, first_order=True)
rho, theta, alpha, drho, dtheta = bb.q_list
drho_old, dtheta_old, dalpha, d2rho, d2theta = bb.u_list

q1 = Symbol('rho')(t)
q2 = Symbol('theta')(t)
q3 = Symbol('alpha')(t)
q4 = q1.diff(t)
q5 = q2.diff(t)
q_names = {rho: q1, theta: q2, alpha: q3, drho: q4, dtheta: q5}
def print_eqns():
    for k, v in eqns.iteritems():
        print str(k) + ': '
        printm(v.subs(q_names))
#print_eqns()

"""
Добавляем новые параметры для электрической части:

La - индуктивность обмотки двигателя
Ra - сопротивление на обмотке двигателя
Kb - constant of the motor’s induced voltage
K2 - torque constant of the motor
U - выходное напряжение (управление),
U0 - значение напряжения на равновесии,
gamma0 - некоторое значение тока, соответствующее положению равновесия системы.
"""
p2 = bb.add_parameters('La Ra Kb K2 U U0 gamma0')
La, Ra, Kb, K2, U, U0, gamma0 = p2

"""
Делаем замену tau = Kg * K2 * (U/Ra - (Kb/Ra)*dtheta),
(U/Ra - (Kb/Ra)*dtheta) - из закона Киркгофа, при La~0.
"""
eqns[d2theta] = eqns[d2theta].subs({tau: -(Kg*K2*Kb/Ra)*dtheta})
#pprint(eqns)

"""
Положение равновесия
"""
q0, u0 = bb.define_equilibrium_point(eqns)
q0[rho] = rho0
q0[theta] = 0
q0[alpha] = 0
q0[drho] = 0
q0[dtheta] = 0

manifold = bb.form_equilibrium_manifold_equations(eqns)
#pprint(manifold)

"""
Численные параметры
"""
p0 = {
    # Ball
    m: 15e-3, r: 9.5e-3, k: 7.0/5, J: 5.41e-7,
    # Beam - steel rod
    l: 400e-3, l1: 160e-3, R: 55e-3, M: 0.3, Jm: 2.66e-2,
    # Electric part
    Ra: 9, La: 0.2e-3, K2: 1.4, Kb: 0.05,
    # Other
    g: 9.8, rho0: 20e-2, Kg: 75
}

"""
Линеаризованные уравнения движения
"""
eqns = linearize(eqns, q0)
#pprint(eqns)
for k in eqns:
    eqns[k] = eqns[k].subs(p0)
#pprint(eqns)

"""
Уравнения возмущенного движения и уравнения первого приближения

peqns = bb.form_perturbed_equations(eqns, manifold)
fa_eqns = bb.form_first_approximation_equations(peqns, q0, params=p0, simplified=False)
#fa_eqns = bb.form_first_approximation_equations(peqns, q0, simplified=False)
#pprint(peqns)
#pprint(fa_eqns)
"""

"""
Матрица коэффициентов уравнений движения системы.
"""

#dx = [x.diff(t) for x in bb.x_list]
#fa_eqns_sorted = [fa_eqns[k] for k in dx]
eqns_sorted = [eqns[u] for u in bb.u_list]
A = bb.create_matrix_of_coeff(eqns_sorted, bb.q_list)
#print A
#exit()

"""
Корни характ. многочлена
"""
eig = A.eigenvals()
for e in eig:
    print e
exit()

B = Matrix([0, 0, 0, 0, 75*10/9])
#pprint(B)

#pprint(ctrb(A,B).tolist())
#if not is_controllable(A, B):
#    print "Pair A,B is not controllable."
    
reg = LQRegulator(A, B)
u = reg.find_control(time=5)
#print u

def f1(x, t):
    # Ограничение на положение, связанное с длинной желоба.
    #print '1: ' + str(x[0])
    #if x[0] < 0:
    #    x[0] = 0.0
    #elif x[0] > 0.4:
    #    x[0] = 0.4
    #print '2: ' + str(x[0])
    dx = A * matrix(x).transpose()
    #print dx[0, 0]
    return mtx2row(dx)

def f2(x, t):
    dx = (A  + B * u) * matrix(x).transpose()
    return mtx2row(dx)

"""
x0 = [0.02, 0, 0, 0, 0]
slv = scipy_odeint(f1, x0, t0=0, t1=50, last=False, h=1e-1)
print slv[:, 1]

t = slv[:, 0]
x1 = slv[:, 1]
x2 = slv[:, 2]
x3 = slv[:, 3]
dx1 = slv[:, 4]
dx2 = slv[:, 5]

plt.figure(1)
plt.plot(t, x1, 'red')
#plt.plot(t, x1, 'red', t, x2, 'green', t, x3, 'blue', t, dx1, 'orange', t, dx2, 'black')
plt.grid(True)
plt.show()
"""

q0 = [0.001, 0.001, 0.001, 0, 0]
slv = scipy_odeint(f2, q0, t0=0, t1=20, last=False, h=1e-1)

t = slv[:, 0]
x1 = slv[:, 1]
x2 = slv[:, 2]
x3 = slv[:, 3]
dx1 = slv[:, 4]
dx2 = slv[:, 5]

plt.figure(1)
plt.plot(t, x1, 'red')
#plt.plot(t, x1, 'red', t, x2, 'green', t, x3, 'blue', t, dx1, 'orange', t, dx2, 'black')
plt.grid(True)
plt.show()



#End