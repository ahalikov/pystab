# coding=windows-1251

"""
Ball&Beam System
"""

from MechanicalFrame import simplify
from MechanicalFrame import normalize
from MechanicalFrame import *

def bb_main():

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
    #printm(eqns[d2theta])
    
    # Добавляю еще одну переменную - силу тока
    current = bb.add_coordinates('q', 1)
    gamma, dgamma, d2gamma = current

    # Добавляю параметры уравнения
    p2 = bb.add_parameters('La Ra Kb K2 U gamma0')
    La, Ra, Kb, K2, U, gamma0 = p2
    """
    La - индуктивность обмотки двигателя
    Ra - сопротивление на обмотке двигателя
    Kb - constant of the motor’s induced voltage
    K2 - torque constant of the motor
    U - выходное напряжение
    gamma0 - некоторое значение тока, соответствующее положению равновесия системы.
    """

    # Делаю замену tau = Kg * K2 * i
    eqns[d2theta] = eqns[d2theta].subs({tau: Kg*K2*gamma})
    #pprint(eqns)
    
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
    manifold = bb.form_equilibrium_manifold_equations(eqns)
    #pprint(manifold)
    gamma0_eqn = solve(manifold[d2theta], gamma0)

    # Численные параметры
    p0 = {
        # Ball
        m: 15e-3, r: 9.5e-3, k: 7.0/5, J: 5.41e-7,
        # Beam - steel rod (стальной стержень)
        l: 400e-3, l1: 160e-3, R: 55e-3, M: 0.5, Jm: 2.66e-2,
        # Electric part
        Ra: 9, La: 0.2e-3, K2: 1e-1,
        # Other
        g: 9.8, rho0: 200e-3, Kg: 75
    }

    p0[gamma0] = gamma0_eqn[0].subs(p0)
    
    # Уравнения возмущенного движения
    peqns = bb.form_perturbed_equations(eqns, manifold)
    #pprint(peqns[d2rho])

    #fa_eqns = bb.form_first_approximation_equations(peqns, q0, simplified=False)
    fa_eqns = bb.form_first_approximation_equations(peqns, q0, params=p0, simplified=False)
    #dx6 = bb.x[dtheta].diff(t)
    #pprint(fa_eqns)

    # Матрица коэффициентов
    A = bb.create_matrix_of_coeff(fa_eqns.values(), bb.x_list)

    # Корни характ. многочлена
    eig = A.eigenvals()
    #pprint(eig)
    

def bb_template():
    """
    Matrix template
    """
    mt = MatrixTemplate()
    q,u,a = mt.define_coordinates(3)
    r,theta,alpha = q
    dr,dtheta,dalpha = u
    T2, T1, T0 = mt.form_ke_matrixes()

    # Making structure as ball and beam system
    for i in range(3):
        for j in range(3):
            T2[i, j] = 0 if i != j else T2[i, j]
        T1[i, 0] = 0
    T2[1,1] = Symbol('a22')(r**2)
    T0 = 0
    #print T2
    mt.set_ke_matrixes(T2, T1, T0)
    #print mt.ke2

    B = mt.form_dhc_matrix(1)
    B[0,0] = 0
    B[0,1] = Symbol('b12')(theta, alpha)
    #print mt.dhc_matrix

    pe = Symbol('U')(r, alpha)
    mt.set_potential_energy(pe)
