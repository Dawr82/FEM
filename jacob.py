from integration import SCHEMA
import itertools

# To mozna zapisac inaczej
# - utworzyc 4 funkcje
# - umiescic te funkcje, w odpowiedniej kolejnosci, w ponizszych tablicach (ilosc funkcji 2 razy mniejsza)

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


# ksi_2w = [[d(p) for d in n_ksi_prim] for p in list(itertools.chain(*zip(SCHEMA[0][1], SCHEMA[0][1])))]
# mu_2w = [[d(p) for d in n_mu_prim] for p in SCHEMA[0][1] + list(reversed(SCHEMA[0][1]))]

# ksi_3w = [[d(p) for d in n_ksi_prim] for p in list(itertools.chain(*zip(SCHEMA[1][1], SCHEMA[1][1]), SCHEMA[1][1]))]
# mu_3w = [[d(p) for d in n_mu_prim] for p in SCHEMA[1][1] + list(reversed(SCHEMA[1][1])) + SCHEMA[1][1]]

# print('\n'.join(['\t'.join([str(n) for n in row]) for row in ksi_2w]))
# print("\n\n")
# print('\n'.join(['\t'.join([str(n) for n in row]) for row in mu_2w]))
# print("\n\n")

# print('\n'.join(['\t'.join([str(n) for n in row]) for row in ksi_3w]))
# print("\n\n")
# print('\n'.join(['\t'.join([str(n) for n in row]) for row in mu_3w]))

# ??? 

# def calculate_matrix(N, derivatives):
#     if N not in (0, 1):
#         raise ValueError("N has to be 0 or 1.")

#     return [[d(p) for d in derivatives] for p in SCHEMA[N][1]]

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


# for matrix in element4_2D.values():
#     print(format_matrix(matrix))

print(*map(format_matrix, element4_2D.values()))




    
        
    






