# coding=utf-8

__author__="Artur"
__date__ ="$01.02.2010 13:40:10$"

from numpy import zeros

def rk45(func, y0, t1, t0=0.0, step=1e-2, whole=True):
        
    y = y0
    t = t0
    h = step

    # Number of equations
    try:
        nequ = len(y0)
    except:
        nequ = 1

    # Number of steps
    steps_num = int(round(t1 - t0) / h) + 1

    # All points of solution
    path = zeros([steps_num, nequ + 1])

    # Function definitions
    int_h = lambda f: h * f
    sum_k = lambda y, k1, k2, k3, k4: y + (k1 + 2*k2 + 2*k3 + k4)/6.0
        
    for i in range(steps_num):

        path[i] = [0.0 for i in range(nequ)]

        if i == 0:
            path[i][0] = t0

        t = t0 + i*h

        # Calculating k1
        f1 = func(y, t)
        k1 = map(int_h, f1)

        # Calculating k2
        f2 = map(lambda yi, ki: func(yi + ki/2.0, t + h/2), y, k1)
        k2 = map(int_h, f2)

        # Calculating k3
        f3 = map(lambda yi, ki: func(yi + ki/2.0, t + h/2), y, k2)
        k3 = map(int_h, f3)

        # Calculating k4
        f4 = map(lambda yi, ki: func(yi + ki, t + h), y, k3)
        k4 = map(int_h, f4)

        # Calculating y
        y = map(sum_k, y, k1, k2, k3, k4)
        path[t] = y

        print t, y
        #print t, h, y, k1, k2, k3, k4
        
        

    return path if whole else y

#End