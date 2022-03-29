import numpy as np

from integration import SCHEMA, N_BC, N_EL
import param


# Shape function derivatives (mu - y axis)
n_ksi_prim = np.array([
    lambda mu : 0.25 * (mu - 1),
    lambda mu : 0.25 * (1 - mu),
    lambda mu : 0.25 * (1 + mu),
    lambda mu : 0.25 * (-1 - mu)
])


# Shape function derivatives (ksi - x axis)
n_mu_prim = np.array([
    lambda ksi : 0.25 * (ksi - 1),
    lambda ksi : 0.25 * (-1 - ksi),
    lambda ksi : 0.25 * (1 + ksi),
    lambda ksi : 0.25 * (1 - ksi)
])


F = np.vectorize(lambda f, x: f(x))

# Shape function derivates' values at integration points. 
SCHEMA_2P_COORDS = np.array(SCHEMA[0][1])
mu_2p = F(n_mu_prim, np.concatenate((SCHEMA_2P_COORDS, SCHEMA_2P_COORDS[::-1]))[np.newaxis].T)
ksi_2p = F(n_ksi_prim, np.concatenate([[n] * 2 for n in SCHEMA[0][1]])[np.newaxis].T)

SCHEMA_3P_COORDS = np.array(SCHEMA[1][1])
mu_3p = F(n_mu_prim, np.concatenate((SCHEMA_3P_COORDS, SCHEMA_3P_COORDS[::-1], SCHEMA_3P_COORDS))[np.newaxis].T)
ksi_3p = F(n_ksi_prim, np.concatenate([[n] * 3 for n in SCHEMA[1][1]])[np.newaxis].T)


# Integration points' weights
SCHEMA_2P_WEIGHTS = np.ones((4, 2))

schema_3p_weights_mu = np.tile(SCHEMA[1][0], 3)[np.newaxis].T
schema_3p_weights_ksi = np.concatenate([[n] * 3 for n in SCHEMA[1][0]])[np.newaxis].T
SCHEMA_3P_WEIGHTS = np.hstack((schema_3p_weights_ksi, schema_3p_weights_mu))


ELEMENT_4_2D = [
    [ksi_2p, mu_2p, SCHEMA_2P_WEIGHTS],
    [ksi_3p, mu_3p, SCHEMA_3P_WEIGHTS],
]


# Jacobian matrix for each integration point inside an element.
def J_matrix_all(nodes, n_schema):
    J_matrices = np.zeros(((n_schema + 2) ** 2, 2, 2))

    if n_schema not in (0, 1):
        return None

    for i, (n_ksi_p, n_mu_p) in enumerate(zip(ELEMENT_4_2D[n_schema][0], ELEMENT_4_2D[n_schema][1])):
        x_ksi = x_mu = y_ksi = y_mu = 0
        for node, n_ksi, n_mu in zip(nodes, n_ksi_p, n_mu_p):
            x_ksi += node.x * n_ksi
            y_ksi += node.y * n_ksi
            x_mu += node.x * n_mu
            y_mu += node.y * n_mu
        J_matrices[i] = np.array([[x_ksi, y_ksi], [x_mu, y_mu]]).round(5)

    return J_matrices


# Jacobian matrix for only one integration point inside an element.
def J_matrix(nodes, n_schema): 
    if n_schema not in (0, 1):
        return None

    x_ksi = x_mu = y_ksi = y_mu = 0
    for node, n_ksi, n_mu in zip(nodes, ELEMENT_4_2D[n_schema][0][0], ELEMENT_4_2D[n_schema][1][0]):
        x_ksi += node.x * n_ksi
        y_ksi += node.y * n_ksi
        x_mu += node.x * n_mu
        y_mu += node.y * n_mu

    return np.array([[x_ksi, y_ksi], [x_mu, y_mu]]).round(5)


def H_matrix(J, n_schema): 
    if n_schema not in (0, 1):
        return None

    H_matrix = np.zeros((4, 4))

    for n_ksi_p, n_mu_p, weights in zip(ELEMENT_4_2D[n_schema][0], ELEMENT_4_2D[n_schema][1], ELEMENT_4_2D[n_schema][2]):
        rj = np.linalg.inv(J)
        dNdx = rj[0, 0] * np.array([n_ksi_p]) + rj[0, 1] * np.array([n_ksi_p])
        dNdy = rj[1, 0] * np.array([n_mu_p]) + rj[1, 1] * np.array([n_mu_p])
        H_matrix += (param.K * (weights[0] * weights[1] * (dNdx.T @ dNdx + dNdy.T @ dNdy)) * np.linalg.det(J)).round(3)
    return H_matrix


def H_matrix_BC(wall_len, wall_id, n_schema):
    if n_schema not in (0, 1):
        return None

    H_BC_matrix = np.zeros((4, 4)) 

    for n_p, weight in zip(N_BC[n_schema][wall_id], SCHEMA[n_schema][0]):
        H_BC_matrix += (weight * n_p * n_p[np.newaxis].T)
    H_BC_matrix *= (param.ALFA * 0.5 * wall_len)
    return H_BC_matrix


def P_vector(wall_len, wall_id, n_schema):
    if n_schema not in (0, 1):
        return None

    P = np.zeros((4, 1))

    for n_p, weight in zip(N_BC[n_schema][wall_id], SCHEMA[n_schema][0]):
        P += weight * n_p.reshape(4, 1)
    P *= (param.T_AMBIENT * param.ALFA * 0.5 * wall_len)
    return P


def C_matrix(J, n_schema):
    if n_schema not in (0, 1):
        return None

    C_mat = np.zeros((4, 4))

    for n_p, weight in zip(N_EL[n_schema], ELEMENT_4_2D[n_schema][2]):
        C_mat +=  (weight[0] * weight[1] * n_p * n_p[np.newaxis].T)
    C_mat*= (param.RO * param.C * np.linalg.det(J))
    return C_mat