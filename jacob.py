from integration import SCHEMA, N_SCHEMA_BC, N_SCHEMA
from constants import *
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
ksi_2w = F(n_ksi_prim, np.ones((4, 4)) * np.concatenate([[n] * 2 for n in SCHEMA[0][1]])[np.newaxis].T)
SCHEMA_2W_WEIGHTS = np.ones((4, 2))


SCHEMA_3W = np.array(SCHEMA[1][1])
mu_3w = F(n_mu_prim, np.ones((9, 4)) * np.concatenate((SCHEMA_3W, SCHEMA_3W[::-1], SCHEMA_3W))[np.newaxis].T)
ksi_3w = F(n_ksi_prim, np.ones((9, 4)) * np.concatenate([[n] * 3 for n in SCHEMA[1][1]])[np.newaxis].T)

schema_3w_weights_mu = np.tile(SCHEMA[1][0], 3)[np.newaxis].T
schema_3w_weights_ksi = np.concatenate([[n] * 3 for n in SCHEMA[1][0]])[np.newaxis].T
SCHEMA_3W_WEIGHTS = np.hstack((schema_3w_weights_ksi, schema_3w_weights_mu))


ELEMENT_4_2D = [
    [ksi_2w, mu_2w, SCHEMA_2W_WEIGHTS],
    [ksi_3w, mu_3w, SCHEMA_3W_WEIGHTS],
]


def J_matrix_all(nodes, N):
    J_matrices = np.zeros(((N + 2) ** 2, 2, 2))

    if N not in (0, 1):
        return None

    for i, (xp, yp) in enumerate(zip(ELEMENT_4_2D[N][0], ELEMENT_4_2D[N][1])):
        x_ksi = x_mu = y_ksi = y_mu = 0
        for node, n_ksi, n_mu in zip(nodes, xp, yp):
            x_ksi += node.x * n_ksi
            y_ksi += node.y * n_ksi
            x_mu += node.x * n_mu
            y_mu += node.y * n_mu
        J_matrices[i] = np.array([[x_ksi, y_ksi], [x_mu, y_mu]]).round(5)

    return J_matrices


def J_matrix(nodes, N): 

    if N not in (0, 1):
        return None

    x_ksi = x_mu = y_ksi = y_mu = 0
    for node, n_ksi, n_mu in zip(nodes, ELEMENT_4_2D[N][0][0], ELEMENT_4_2D[N][1][0]):
        x_ksi += node.x * n_ksi
        y_ksi += node.y * n_ksi
        x_mu += node.x * n_mu
        y_mu += node.y * n_mu

    return np.array([[x_ksi, y_ksi], [x_mu, y_mu]]).round(5)


def H_matrix(J, N): 

    if N not in (0, 1):
        return None

    H_matrix = np.zeros((4, 4))

    for nx, ny, weights in zip(ELEMENT_4_2D[N][0], ELEMENT_4_2D[N][1], ELEMENT_4_2D[N][2]):
        rj = np.linalg.inv(J)
        dNdx = rj[0, 0] * np.array([nx]) + rj[0, 1] * np.array([nx])
        dNdy = rj[1, 0] * np.array([ny]) + rj[1, 1] * np.array([ny])
        H_matrix += (K * (weights[0] * weights[1] * (dNdx.T @ dNdx + dNdy.T @ dNdy)) * np.linalg.det(J)).round(3)
    return H_matrix


def H_matrix_BC(L, wall_id, N):
    if N not in (0, 1):
        return None

    H_BC = np.zeros((4, 4)) 

    for point, weight in zip(N_SCHEMA_BC[N][wall_id], SCHEMA[N][0]):
        H_BC += (weight * point * point[np.newaxis].T)
    H_BC *= (ALFA * 0.5 * L)
    return H_BC


def P_vector(L, wall_id, N):
    if N not in (0, 1):
        return None

    P = np.zeros((4, 1))

    for point, weight in zip(N_SCHEMA_BC[N][wall_id], SCHEMA[N][0]):
        P += weight * point.reshape(4, 1)
    P *= (T_AMBIENT * ALFA * 0.5 * L)
    return P


def C_matrix(J, N):
    if N not in (0, 1):
        return None

    C_mat = np.zeros((4, 4))

    for N, weight in zip(N_SCHEMA[N], ELEMENT_4_2D[N][2]):
        C_mat +=  weight[0] * weight[1] * N * N[np.newaxis].T
    C_mat*= (RO * C * np.linalg.det(J))
    return C_mat



    
        
    






