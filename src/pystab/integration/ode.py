# coding=utf-8

__author__="Artur"
__date__ ="$01.02.2010 13:40:10$"

from numpy import zeros, linspace

def step_rkf45(derivs, x, t, h):
    """
    Шаг метода Рунге-Кутты 4 порядка.
    """
    int_h = lambda f: h * f
    sum_k = lambda x, k1, k2, k3, k4: x + (k1 + 2*k2 + 2*k3 + k4)/6.0

	# Calculating k1, k2, k3, k4 and k5
    try:
        k1 = map(int_h, derivs(x, t))
        k2 = map(int_h, derivs(map(lambda xi, ki: xi + ki/2.0, x, k1), t + h/2))
        k3 = map(int_h, derivs(map(lambda xi, ki: xi + ki/2.0, x, k2), t + h/2))
        k4 = map(int_h, derivs(map(lambda xi, ki: xi + ki, x, k3), t + h))
    except:
        print x

    # Calculating x
    x = map(sum_k, x, k1, k2, k3, k4)
    
    return x

def odeint(derivs, x0, t1, t0=0.0, h=1e-3, step=step_rkf45, last=True):
    """
    Численное интегрирование ОДУ.
    """
    try:
        nequ = len(x0)
    except:
        nequ = 1

    x = x0
    n = int(round(t1 - t0) / h)

    slv = zeros([n + 1, nequ + 1])
    slv[0, 0] = t0
    slv[0, 1:] = x

    for i in range(1, n + 1):
        t = t0 + i*h
        x = step(derivs, x, t, h)
        slv[i, 0] = t
        slv[i, 1:] = x

    return x if last else slv

def scipy_odeint(derivs, x0, t1, t0=0.0, h=1e-3, last=True):
    """
    Обертка для функции odeint из библиотеки LAPACK(scipy) для числ. интегрирования ОДУ.
    Работает быстрее и безглючнее, чем самописная odeint выше, но для использования нужно
    установить библиотеку scipy (scipy.org).
    """
    from scipy.integrate import odeint
    n = int(round((t1 - t0)/h))
    t = linspace(0, t1, n)
    x = odeint(derivs, x0, t)

    # If only last value of solution is needed
    if last:
        return x[-1]

    # Let's make complete solution array: time + trac
    slv = zeros([n, x.shape[1] + 1])
    for i in range(n):
        slv[i, 0] = t[i]
        slv[i, 1:] = x[i, :]
        
    return slv

"""
For quick tests
"""
def main():
    # x' = x
    f = lambda x, t: x
    s = odeint(f, [1], 1, h=1e-1)
    print s

    print scipy_odeint(f, [1], 1)

if __name__ == '__main__':
    main()

#End