from math import sqrt
from inspect import signature
import numpy as np
import itertools


class IntegrationError(Exception):
    pass


a = 1/sqrt(3)
b = sqrt(0.6)


# Kwadratura Gaussa (2 i 3 punktowe schematy)
SCHEMA = [
    [[1, 1],[-a, a]],                                                       
    [[5/9, 8/9, 5/9],[-b, 0, b]],                                             
]


# Calkowanie numeryczne Gaussa (przekazac funkcję (1 lub 2 argumenty) oraz schemat calkowania)
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

# Funkcje ksztaltu dla elementu 2D 4 wezlowego
N = [
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 - mu),
    lambda ksi, mu : 0.25 * (1 + ksi) * (1 + mu),
    lambda ksi, mu : 0.25 * (1 - ksi) * (1 + mu),
]

# Wspolrzedne punktow na scianach elementu w 2-punktowym schemacie calkowania
SCHEMA_2P_BC = np.array([
    [(-a, -1), (a, -1)],      
    [(1, -a), (1, a)],    
    [(a, 1), (-a, 1)],    
    [(-1, a), (-1, -a)], 
])

# Wspolrzedne punktow na scianach elementu w 3-punktowym schemacie calkowania
SCHEMA_3P_BC = np.array([
    [(-b, -1), (0, -1), (b, -1)],
    [(1, -b), (1, 0), (1, b)],
    [(b, 1), (0, 1), (-b, 1)],
    [(-1, b), (-1, 0), (-1, -b)],
])


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


# Wektory wartosci wszystkich funkcji ksztaltu w kazdym punkcie na kazdej scianie (2 i 3 punktowy schemat calkowania)
# Do liczenia macierzy Hbc oraz wektora P.
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


# Wartosci funkcji ksztaltu w punktach calkowania 
# Do liczenia macierzy C.
N_EL = [
    n_schema_el_init(SCHEMA[0][1]),
    n_schema_el_init(SCHEMA[1][1]),
]

