# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

"""
Ball&Beam System
"""

from pystab.mechanics import expand
from sympy import latex

from pystab.stability import *
from pystab.mechanics import *

bb = MechanicalFrame(name='Ball on the Beam')

"""
Координаты
rho - положение шарика на желобе, отсчитываемое от левого неподвижного
конца желоба;
theta - угловая скорость колеса двигателя;
alpha - угол на который поднимается желоб.
"""
q, u, a = bb.add_coordinates('q', 3)
rho, theta, alpha = q
drho, dtheta, dalpha = u
d2rho, d2theta, d2alpha = a

# Параметры системы
p1 = bb.add_parameters('m r J1 M l J2 l1 R J3 g k Kg rho0 theta0 alpha0 b0 tau')
m, r, J1, M, l, J2, l1, R, J3, g, k, Kg, rho0, theta0, alpha0, b0, tau = p1
"""
m, r, J1 - масса, радиус и момент инерции шарика
M, l, J2 - масса, длина и момент инерции желоба
l1 - длина рычага, соединяющего колесо двигателя с желобом
R, J3 - радиус и момент инерции колеса
g - гравитационная постоянная
k = (5 + 2*(R/R')^2)/5, где R' - малый радиус качения шарика.
Kg - передаточное число
b0 - коэффициент трения (диссипация)
rho0 - исследуемое положение шарика
tau - управляющий момент на двигателе (tau = Kg * K2 * i)
где K2 - это постоянная момента двигателя, i - ток.
"""

# Обобщенные силы
Q = bb.add_joint_forces({rho: 0, theta: -b0*dtheta + tau, alpha: 0})

# Кинетическая энергия
T = 1.0/2*((J1 + m)*drho**2 + (J2 + m*(rho**2))*dalpha**2 + J3*dtheta**2)
# Потенциальная энергия
P = -g*sin(alpha)*(m*rho + 1.0/2*M*l)
# Функция Лагранжа
L =  T - P
#pprint(L)
bb.set_lagrangian(L)

# Уравнение связи
dhc_eqn_full = [l**2*dalpha*sin(alpha) + R**2*dtheta*sin(theta)+ \
    l*R*(dalpha - dtheta)*sin(alpha - theta) - l*R*(dtheta*sin(theta)+ \
    dalpha*sin(alpha)) + l*l1*dalpha*cos(alpha) - l1*R*dtheta*cos(theta)]

# Уравнение связи при alpha и theta ~ 0
dhc_eqn_simple = [l*dalpha - R*dtheta]
dhc_eqn = dhc_eqn_full
#dhc_eqn = dhc_eqn_simple

# Уравнение связи, разрешенное относительно зависимой скорости
dhc_eqn = solve(dhc_eqn, [dalpha])
dhc_eqn = collect(dhc_eqn[dalpha], dtheta)
bb.form_constraints_matrix([dhc_eqn], [dalpha])

# Уравнения Шульгина для r, theta, alpha
eqns = bb.form_shulgin_equations(normalized=True)
#pprint(eqns[d2rho])
#exit()

# Добавляю еще одну переменную - силу тока
current = bb.add_coordinates('q', 1)
gamma, dgamma, d2gamma = current

# Добавляю параметры уравнения
p2 = bb.add_parameters('La Ra Kb K2 U U0 gamma0')
La, Ra, Kb, K2, U, U0, gamma0 = p2
"""
La - индуктивность обмотки двигателя
Ra - сопротивление на обмотке двигателя
Kb - constant of the motor’s induced voltage
K2 - torque constant of the motor
U - выходное напряжение (управление),
U0 - значение напряжения на равновесии,
gamma0 - некоторое значение тока, соответствующее положению равновесия системы.
"""

# Делаю замену tau = Kg * K2 * i
eqns[d2theta] = eqns[d2theta].subs({tau: Kg*K2*gamma})
#pprint(eqns[d2theta])

# Добавляю уравнение для тока
current_eqn = La*dgamma + Ra*gamma + Kb*dtheta - U
eqns[dgamma] = U/La - (Ra/La)*gamma - (Kb/La)*dtheta
#pprint(current_eqn)

# Положение равновесия
q0, u0 = bb.define_equilibrium_point(eqns)
q0[rho] = rho0
q0[theta] = 0
q0[alpha] = 0
q0[gamma] = gamma0
#print q0, u0
manifold = bb.form_equilibrium_manifold_equations(eqns, {U: U0})
#print 'Equilibrium point: ', manifold
gamma0_eqn = solve(manifold[dgamma], gamma0)

# Численные параметры
p0 = {
    # Шарик
    m: 15e-3, r: 9.5e-3, k: 7.0/5, J1: 5.41e-7,
    # Желоб и рычаг
    l: 400e-3, l1: 160e-3, M: 0.3, R: 55e-3, J2: 1.6e-2,
    # Колесо двигателя
    J3: 0.11e-2,
    # Электрическая часть
    Ra: 9, La: 0.2e-3, K2: 1.4, Kb: 1,
    # Другие параметры
    g: 9.8, rho0: 200e-3, Kg: 75.0,
    #
    U0: 0,
    b0: 0.1
}

p0[gamma0] = gamma0_eqn[0].subs(p0)

# Уравнения возмущенного движения
peqns = bb.form_perturbed_equations(eqns, manifold)

#fa_eqns = bb.form_approximated_equations(peqns, q0, simplified=False)
fa_eqns = bb.form_approximated_equations(peqns, q0, params=p0, simplified=False)
#dx6 = bb.x[dtheta].diff(t)
#printm(fa_eqns)
#exit()

x_0 = [(x, 0) for x in bb.x.values()]

#eqn = peqns[d2theta]
#tmp = (pdiff(eqn, x3).subs({x1:0, x3:0, x5:0, x6:0}))
#pprint(tmp)

# Матрица коэффициентов
dx =  [x.diff(t) for x in bb.x_list]
fa_eqns_sorted = [fa_eqns[k] for k in dx]
A = bb.create_matrix_of_coeff(fa_eqns_sorted, bb.x_list)
#A[5, 0] = Symbol('A61')
#A[5, 3] = Symbol('A64')
pprint(A)
#print A.tolist()

# Корни характ. многочлена
eig = A.eigenvals()
pprint(eig)

#B = Matrix([0, 0, 0, 1/0.2e-3, 0, 0])
#B = Matrix([0, 0, 0, 1/La, 0, 0])
#B = Matrix([0, 0, 0, 1, 0, 0])
#pprint(B)

#C = sym_ctrb(A, B)
#print C.det()

#reg = LQRegulator(A, B)
#u = reg.find_control(time=10)
#print u