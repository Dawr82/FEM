from math import sqrt
from inspect import signature
import itertools

import numpy as np


class IntegrationError(Exception):
    """Raised where integration-related errors occur."""
    pass


a = 1/sqrt(3)
b = sqrt(0.6)


# Gaussian Quadrature (2 and 3 point schemas)
SCHEMA = [
    [[1, 1],[-a, a]],                                                       
    [[5/9, 8/9, 5/9],[-b, 0, b]],                                             
]


# Gaussian integration (for 2 and 3 point schemas)
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


#=======================================================================================

# Shape functions for 2D 4-node element.
N = [
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 + mu),
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 + mu),
]

# Point coordinates at element's walls - 2-point Gaussian integration
SCHEMA_2P_BC = np.array([
    [(-a, -1), (a, -1)],      
    [(1, -a), (1, a)],    
    [(a, 1), (-a, 1)],    
    [(-1, a), (-1, -a)], 
])

# Point coordinates at element's walls - 3-point Gaussian integration
SCHEMA_3P_BC = np.array([
    [(-b, -1), (0, -1), (b, -1)],
    [(1, -b), (1, 0), (1, b)],
    [(b, 1), (0, 1), (-b, 1)],
    [(-1, b), (-1, 0), (-1, -b)],
])


# Helper function used to set up shape function vectors.
def n_schema_bc_init(schema_bc):
    n_schema_bc = np.zeros((4, len(schema_bc[0]), 4))
    for i, wall in enumerate(schema_bc):
        for j, point in enumerate(wall):
            if i == len(SCHEMA_2P_BC) - 1:
                n_schema_bc[i, j, i] = N[i](*point)
                n_schema_bc[i, j, 0] = N[0](*point)
            else:
                n_schema_bc[i, j, i] = N[i](*point)
                n_schema_bc[i, j, i + 1] = N[i + 1](*point)
    return n_schema_bc


# Shape function values at each integration point of all element's walls.
N_BC = [
    n_schema_bc_init(SCHEMA_2P_BC),
    n_schema_bc_init(SCHEMA_3P_BC),
]


def n_schema_el_init(schema_el):
    n_schema_el = np.zeros((len(schema_el) ** 2 ,4))
    for i, (ksi, mu) in enumerate(list(itertools.product(schema_el, schema_el))):
        for j, n in enumerate(N):
            n_schema_el[i][j] = n(ksi, mu)
    return n_schema_el


# Shape function values at integration points of an element (inside the element).
N_EL = [
    n_schema_el_init(SCHEMA[0][1]),
    n_schema_el_init(SCHEMA[1][1]),
]