from math import sqrt
from inspect import signature
import numpy as np
import itertools


class IntegrationError(Exception):
    pass


SCHEMA = [
    [[1, 1],[-1/sqrt(3), 1/sqrt(3)]],                                                       
    [[5/9, 8/9, 5/9],[-sqrt(3/5), 0, sqrt(3/5)]],                                             
]


def gauss_integration(f, N):
    n_vars = len(signature(f).parameters)
    result = 0

    try:
        if n_vars == 1:
            for i in range(N + 2):
                result += f(SCHEMA[N][1][i]) * SCHEMA[N][0][i]
        elif n_vars == 2:
            for i in range(N + 2):
                for j in range(N + 2):
                    result += f(SCHEMA[N][1][i], SCHEMA[N][1][j]) * SCHEMA[N][0][i] * SCHEMA[N][0][j]
        else:
            raise IntegrationError("Function has to have 1 or 2 variables.")
    except IndexError:
        print("SCHEMA doesn't support that N.")

    return result


a = 1/sqrt(3)
b = sqrt(3/5)

SCHEMA_BC_2W = [
    [(-a, -1), (a, -1)],      
    [(1, -a), (1, a)],    
    [(a, 1), (-a, 1)],    
    [(-1, a), (-1, -a)], 
]

SCHEMA_BC_3W = [
    [(-b, -1), (0, -1), (b, -1)],
    [(1, -b), (1, 0), (1, b)],
    [(b, 1), (0, 1), (-b, 1)],
    [(-1, b), (-1, 0), (-1, -b)],
]


N = [
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 + mu),
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 + mu),
]

def n_schema_bc_init(schema_bc):
    n_schema_bc = np.zeros((4, len(schema_bc[0]), 4))
    for i, wall in enumerate(schema_bc):
        for j, point in enumerate(wall):
            if i == len(SCHEMA_BC_2W) - 1:
                n_schema_bc[i, j, i] = N[i](*point)
                n_schema_bc[i, j, 0] = N[0](*point)
            else:
                n_schema_bc[i, j, i] = N[i](*point)
                n_schema_bc[i, j, i + 1] = N[i + 1](*point)
    return n_schema_bc

N_SCHEMA_BC = [
    n_schema_bc_init(SCHEMA_BC_2W),
    n_schema_bc_init(SCHEMA_BC_3W),
]


N_SCHEMA_2W = np.zeros((4, 4))
N_SCHEMA_3W = np.zeros((9, 4))

for i, (ksi, mu) in enumerate(list(itertools.product(SCHEMA[0][1], SCHEMA[0][1]))):
    for j, n in enumerate(N):
        N_SCHEMA_2W[i][j] = n(ksi, mu)


for i, (ksi, mu) in enumerate(list(itertools.product(SCHEMA[1][1], SCHEMA[1][1]))):
    for j, n in enumerate(N):
        N_SCHEMA_3W[i][j] = n(ksi, mu)


N_SCHEMA = [
    N_SCHEMA_2W,
    N_SCHEMA_3W,
]