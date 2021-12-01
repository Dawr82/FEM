from math import sqrt
import itertools

N = 100

def rec_integration(f, a, b, n):

    sum = 0.0
    dx = (b - a) / n

    for i in range(n):
        x0 = a + i * dx
        x1 = x0 + dx
        sum += f((x0 + x1) / 2) * dx

    return sum


def trap_integration(f, a, b, n):

    sum = 0.0
    dx = (b - a) / n

    for i in range(n):
        x0 = a + i * dx
        x1 = x0 + dx
        sum += 0.5 * dx * (f(x0) + f(x1))

    return 
    

def gauss_integration_1D(f, n):
    # n = 1 -> 2 punkty calkowania (oraz 2 wagi z nimi zwiazane)
    # n = 2 -> 3 punkty calkowania (oraz 3 wagi z nimi zwiazane)

    if n == 1:
        w = (1.0, 1.0)
        xc = (-1 / sqrt(3), 1 / sqrt(3))
    elif n == 2:
        w = (5/9, 8/9, 5/9)
        xc = (-sqrt(3/5), 0, sqrt(3/5))
    else:
        return None
        
    sum = 0.0
    for w_i, xc_i in zip(w, xc):
        sum += f(xc_i) * w_i

    return sum


def gauss_integration_2D(f, n):
    # n = 1 lub n = 2 else error

    if n == 1:
        w = (1.0, 1.0)
        xc = (-1 / sqrt(3), 1 / sqrt(3))
    elif n == 2:
        w = (5/9, 8/9, 5/9)
        xc = (-sqrt(3/5), 0, sqrt(3/5))
    else:
        return None

    print(list(itertools.product(w, w)))
    print(list(itertools.product(xc, xc)))


    #for w_i, xc_i in zip(itertools.product(w, w), itertools.product(xc, xc)):




        


f = lambda x : x ** 2

print(rec_integration(f, -3.0, 3.0, N))
print(trap_integration(f, -3.0, 3.0, N))

gauss_integration_2D(f, 2)



