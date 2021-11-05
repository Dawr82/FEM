from integration import SCHEMA
import itertools

# TODO: try to replace these two list with one (both contain the same functions but in different order)

n_ksi_prim = [
    lambda mu : 0.25 * (mu - 1),
    lambda mu : 0.25 * (1 - mu),
    lambda mu : 0.25 * (1 + mu),
    lambda mu : 0.25 * (-1 - mu),
]


n_mu_prim = [
    lambda ksi : 0.25 * (ksi - 1),
    lambda ksi : 0.25 * (-1 - ksi),
    lambda ksi : 0.25 * (1 + ksi),
    lambda ksi : 0.25 * (1 - ksi),
]


mu_2w = [[d(p) for d in n_ksi_prim] for p in SCHEMA[0][1]]
ksi_2w = [[d(p) for d in n_mu_prim] for p in SCHEMA[0][1]]

mu_3w = [[d(p) for d in n_ksi_prim] for p in SCHEMA[1][1]]
ksi_3w = [[d(p) for d in n_mu_prim] for p in SCHEMA[1][1]]


element4_2D = {
    "ksi2w" : list(itertools.chain(*zip(mu_2w, mu_2w))),
    "mu2w" : ksi_2w + ksi_2w[::-1],
    "ksi3w" : list(itertools.chain(*zip(mu_3w, mu_3w, mu_3w))),
    "mu3w" : ksi_3w + ksi_3w[::-1] + ksi_3w,
}

def format_matrix(matrix):
    return '\n' + '\n'.join(['\t'.join([str(n) for n in row]) for row in matrix]) + '\n'


#print(*map(format_matrix, element4_2D.values()))




    
        
    






