from math import sqrt
from inspect import signature
import numpy as np
import itertools


class IntegrationError(Exception):
    pass

SCHEMA = [
    [[1, 1],[-1/sqrt(3), 1/sqrt(3)]],                                                            # 2 w
    [[5/9, 8/9, 5/9],[-sqrt(3/5), 0, sqrt(3/5)]],                                                # 3 w
    [[0.347855, 0.652145, 0.652145,  0.347855],[-0.861136, -0.339981, 0.339981, 0.861136]],      # 4 w
]


SCHEMA_BC_2W = [
    [(-1/sqrt(3), -1), (1/sqrt(3), -1)],  # 1 ściana      
    [(1, -1/sqrt(3)), (1, 1/sqrt(3))],    # 2 ściana
    [(1/sqrt(3), 1), (-1/sqrt(3), 1)],    # 3 ściana
    [(-1, 1/sqrt(3)), (-1, -1/sqrt(3))],  # 4 ściana
]

N = [
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 + mu),
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 + mu),
]


N_SCHEMA_2W = np.zeros((4, 4))
N_SCHEMA_3W = np.zeros((9, 4))

WEIGHT_SCHEMA_2W = np.array(list(itertools.product(SCHEMA[0][0], SCHEMA[0][0])))
WEIGHT_SCHEMA_3W = np.array(list(itertools.product(SCHEMA[1][0], SCHEMA[1][0])))

for i, (ksi, mu) in enumerate(list(itertools.product(SCHEMA[0][1], SCHEMA[0][1]))):
    for j, n in enumerate(N):
        N_SCHEMA_2W[i][j] = n(ksi, mu)


for i, (ksi, mu) in enumerate(list(itertools.product(SCHEMA[1][1], SCHEMA[1][1]))):
    for j, n in enumerate(N):
        N_SCHEMA_3W[i][j] = n(ksi, mu)


N_SCHEMA_BC_2W = np.zeros((4, 2, 4))

for i, wall in enumerate(SCHEMA_BC_2W):
    for j, point in enumerate(wall):
        if i == len(SCHEMA_BC_2W) - 1:
            N_SCHEMA_BC_2W[i, j, i] = N[i](*point)
            N_SCHEMA_BC_2W[i, j, 0] = N[0](*point)
        else:
            N_SCHEMA_BC_2W[i, j, i] = N[i](*point)
            N_SCHEMA_BC_2W[i, j, i + 1] = N[i + 1](*point)


def integral_gauss(f, N):
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



