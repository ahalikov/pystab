# coding=utf-8

__author__="Artur"
__date__ ="$26.01.2010 16:07:56$"

"""
Ball&Beam System
"""

from pystab.stability import *
from pystab.mechanics import *

bb = MechanicalFrame(name='Ball on the Beam')

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

# Параметры системы
p1 = bb.add_parameters('m r J M l l1 R Jm g k Kg rho0 tau')
m, r, J, M, l, l1, R, Jm, g, k, Kg, rho0, tau = p1
"""
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

# Обобщенные силы
Q = bb.add_joint_forces({r: 0, theta: tau, alpha: 0})

# Кинетическая энергия
T = 1.0/2*(k*m*drho**2 + (J + m*(rho**2))*dalpha**2 + Jm*dtheta**2)
# Потенциальная энергия
P = g*sin(alpha)*(m*rho + 1.0/2*M*l)
# Функция Лагранжа
L =  T - P
#pprint(L)
bb.set_lagrangian(L)

# Уравнение связи
dhc_eqn_full = [l**2*dalpha*sin(alpha) + R**2*dtheta*sin(theta)+ \
    l*R*(dalpha - dtheta)*sin(alpha - theta) - l*R*(dtheta*sin(theta)+ \
    dalpha*sin(alpha)) + l*l1*dalpha*cos(alpha) - l1*R*dtheta*cos(theta)]

# Уравнение связи когда alpha и theta ~ 0
dhc_eqn_simple = [l*dalpha - R*dtheta]

dhc_eqn = dhc_eqn_full

# Уравнение связи, разрешенное относительно зависимой скорости
dhc_eqn = solve(dhc_eqn, [dalpha])
dhc_eqn = collect(dhc_eqn[dalpha], dtheta)
bb.form_constraints_matrix([dhc_eqn], [dalpha])

# Уравнения Шульгина для r, theta, alpha
eqns = bb.form_shulgins_equations(normalized=True, expanded=False)
#pprint(eqns[d2theta])

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
#print "Equilibrium point:"
#pprint(manifold)
gamma0_eqn = solve(manifold[d2theta], gamma0)
#pprint(gamma0_eqn)

# Численные параметры
p0 = {
    # Ball
    m: 15e-3, r: 9.5e-3, k: 7.0/5, J: 5.41e-7,
    # Beam - steel rod
    l: 400e-3, l1: 160e-3, R: 55e-3, M: 0.5, Jm: 2.66e-2,
    # Electric part
    Ra: 9, La: 0.2e-3, K2: 1.4, Kb: 1e-1,
    # Other
    g: 9.8, rho0: 200e-3, Kg: 75.0
}

p0[gamma0] = gamma0_eqn[0].subs(p0)
#print p0[gamma0]

# Уравнения возмущенного движения
peqns = bb.form_perturbed_equations(eqns, manifold)
#pprint(peqns[dgamma])

#fa_eqns = bb.form_first_approximation_equations(peqns, q0, simplified=False)
#fa_eqns = bb.form_first_approximation_equations(peqns, q0, params=p0, simplified=False)
#dx6 = bb.x[dtheta].diff(t)
#pprint(fa_eqns)

x1 = bb.x[rho]
x2 = bb.x[theta]
x3 = bb.x[alpha]
x5 = bb.x[drho]
x6 = bb.x[dtheta]

printm(bb.reduced_lagrangian)
tmp1 = bb.dhc_matrix[0, 1]*pdiff(bb.reduced_lagrangian, alpha)
for q in bb.q_list:
    u = q.diff(t)
    if u in bb.x:
        tmp1 = tmp1.subs(u, bb.x.get(u) + bb.u0.get(u))
    tmp1 = tmp1.subs(q, bb.x.get(q) + bb.q0.get(q))

tmp1 = pdiff(tmp1, x2)

#printm(tmp1)

x_0 = [(x, 0) for x in bb.x.values()]

eqn = peqns[d2theta]
#tmp = (pdiff(eqn, x3).subs({x1:0, x3:0, x5:0, x6:0}))
#pprint(tmp)



# Матрица коэффициентов
#dx =  [x.diff(t) for x in bb.x_list]
#fa_eqns_sorted = [fa_eqns[k] for k in dx]
#A = bb.create_matrix_of_coeff(fa_eqns_sorted, bb.x_list)
#pprint(A)
#print A.tolist()

# Корни характ. многочлена
#eig = A.eigenvals()
#pprint(eig)

#B = Matrix([0, 0, 0, 1/0.2e-3, 0, 0])
#B = Matrix([0, 0, 0, 1/La, 0, 0])
#pprint(B)

#C = ctrb(A, B)
#pprint(C)

#reg = LQRegulator(A, B)
#u = reg.find_control(time=5)
#print u