from integration import SCHEMA
import itertools
import numpy as np


n_ksi_prim = np.array([
    lambda mu : 0.25 * (mu - 1),
    lambda mu : 0.25 * (1 - mu),
    lambda mu : 0.25 * (1 + mu),
    lambda mu : 0.25 * (-1 - mu)
])


n_mu_prim = np.array([
    lambda ksi : 0.25 * (ksi - 1),
    lambda ksi : 0.25 * (-1 - ksi),
    lambda ksi : 0.25 * (1 + ksi),
    lambda ksi : 0.25 * (1 - ksi)
])

F = np.vectorize(lambda f, x: f(x))

SCHEMA_2W = np.array(SCHEMA[0][1])
mu_2w = F(n_mu_prim, np.ones((4, 4)) * np.concatenate((SCHEMA_2W, SCHEMA_2W[::-1]))[np.newaxis].T)
ksi_2w = F(n_ksi_prim, np.ones((4, 4)) * np.tile(SCHEMA_2W,  2)[np.newaxis].T)

SCHEMA_3W = np.array(SCHEMA[1][1])
mu_3w = F(n_mu_prim, np.ones((9, 4)) * np.concatenate((SCHEMA_3W, SCHEMA_3W[::-1], SCHEMA_3W))[np.newaxis].T)
ksi_3w = F(n_ksi_prim, np.ones((9, 4)) * np.tile(SCHEMA_3W,  3)[np.newaxis].T)


# ksi_2w = [[d(p) for d in n_ksi_prim] for p in SCHEMA[0][1]]
# mu_2w = [[d(p) for d in n_mu_prim] for p in SCHEMA[0][1]]

# ksi_3w = [[d(p) for d in n_ksi_prim] for p in SCHEMA[1][1]]
# mu_3w = [[d(p) for d in n_mu_prim] for p in SCHEMA[1][1]]


element4_2D = {
    # "ksi2w" : list(itertools.chain(*zip(ksi_2w, ksi_2w))),
    # "mu2w" : mu_2w + mu_2w[::-1],
    # "ksi3w" : list(itertools.chain(*zip(ksi_3w, ksi_3w, ksi_3w))),
    # "mu3w" : mu_3w + mu_3w[::-1] + mu_3w,
    "ksi2w" : ksi_2w,
    "mu2w" : mu_2w,
    "ksi3w" : ksi_3w,
    "mu3w" : mu_2w,
}

def format_matrix(matrix):
    return '\n' + '\n'.join(['\t'.join([str(n) for n in row]) for row in matrix]) + '\n'


#print(*map(format_matrix, element4_2D.values()))




    
        
    






