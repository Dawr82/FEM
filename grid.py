from jacob import element4_2D
from integration import *
from globals import *
import numpy as np
import itertools

np.set_printoptions(linewidth=200)


def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)


def J_matrix_all(nodes):
    jacobians = np.zeros((4, 2, 2))

    for i, (xp, yp) in enumerate(zip(element4_2D.get('ksi2w'), element4_2D.get('mu2w'))):
        x_ksi = x_mu = y_ksi = y_mu = 0
        for node, n_ksi, n_mu in zip(nodes, xp, yp):
            x_ksi += node.x * n_ksi
            y_ksi += node.y * n_ksi
            x_mu += node.x * n_mu
            y_mu += node.y * n_mu
        jacobians[i] = np.array([[x_ksi, y_ksi], [x_mu, y_mu]]).round(4)

    return jacobians


def J_matrix(nodes):
    x_ksi = x_mu = y_ksi = y_mu = 0
    for node, n_ksi, n_mu in zip(nodes, element4_2D.get('ksi2w')[0], element4_2D.get('mu2w')[0]):
        x_ksi += node.x * n_ksi
        y_ksi += node.y * n_ksi
        x_mu += node.x * n_mu
        y_mu += node.y * n_mu
    #return np.array([x_ksi, y_ksi, x_mu, y_mu]).round(7)
    return [[[round(x_ksi, 6), round(y_ksi, 6)], [round(x_mu, 6), round(y_mu, 6)]]] * 4


def H_matrix(J): 
    H_matrix = np.zeros((4, 4))
    for nx, ny, j in zip(element4_2D.get('ksi2w'), element4_2D.get('mu2w'), J):
        rj = np.linalg.inv(j)
        dNdx = rj[0, 0] * np.array([nx]) + rj[0, 1] * np.array([nx])
        dNdy = rj[1, 0] * np.array([ny]) + rj[1, 1] * np.array([ny])
        H_matrix += (K * (dNdx.T @ dNdx + dNdy.T @ dNdy) * np.linalg.det(j)).round(3)
    return H_matrix


def H_matrix_BC(J, wall_id):
    H_BC = np.zeros((4, 4)) 
    for point in N_SCHEMA_BC_2W[wall_id]:
        H_BC += (point * point[np.newaxis].T)
    H_BC *= (ALFA * np.linalg.det(J))
    return H_BC


def P_vector(J, wall_id):
    P = np.zeros((4, 1))
    for point in N_SCHEMA_BC_2W[wall_id]:
        print(point, point.shape)
        P += point.reshape(4, 1)
    P *= (T_AMBIENT * ALFA * np.linalg.det(J))
    return P


def C_matrix(J):
    C_mat = np.zeros((4, 4))
    for N in trunc(N_SCHEMA_2W, 4):
        C_mat += N * N[np.newaxis].T
    C_mat*= (RO * C * np.linalg.det(J))
    return C_mat
    

class Node:

    def __init__(self, x, y, t0, bc=0):
        self.x = x
        self.y = y
        self.bc = bc
        self.t0 = t0

    def __repr__(self):
        return f"[{self.x: .3f}, {self.y: .3f}, BC = {self.bc}]" 


class Element:

    def __init__(self, ids, jacobians=None, H_matrix=None, H_BC_matrices=None, H_BC_matrix=None, P_vectors=None, P_vector=None, C_matrix=None):
        self.ids = ids 
        self.jacobians = jacobians
        self.H_matrix = H_matrix
        self.H_BC_matrices = H_BC_matrices
        self.H_BC_matrix = H_BC_matrix
        self.P_vectors = P_vectors
        self.P_vector = P_vector
        self.C_matrix = C_matrix


    def __repr__(self):
        return 30 * "=" + f"\n\n{self.ids}\n\nJ matrix:\n\n{self.jacobians}\n\nH matrix:\n\n{self.H_matrix}\n\n \
H_BC_matrix\n\n{self.H_BC_matrix}\n\P_vector\n\n{self.P_vector}\n\nC_matrix\n\n{self.C_matrix}\n\n"


    def __len__(self):
        return len(self.ids)

  
class Grid:

    def __init__(self, height, breadth, num_nodes_height, num_nodes_breadth, temp_start):
        self.h = height
        self.b = breadth
        self.n_b = num_nodes_breadth
        self.n_h = num_nodes_height
        self.n_n = self.n_h * self.n_b
        self.n_e = (self.n_h - 1) * (self.n_b - 1)
        self.t0 = temp_start
        self.nodes = self.create_nodes()
        self.elements = self.create_elements()
  
        self.calculate_jacobians()
        self.calculate_H_matrices()
        self.apply_boundary_conditions()
        self.update_H_matrices()
        self.calculate_C_matrices()

        self.H_global = self.aggregate_H()
        self.P_global = self.aggregate_P()
        self.C_global = self.aggregate_C()


    def set_boundary_condition(self, nodes):
        for node in nodes:
            if any((node.x == 0, node.y == 0, node.x == self.b, node.y == self.h)):
                node.bc = 1


    def create_nodes(self):
        nodes = [Node(j * (self.b / (self.n_b - 1)), i * (self.h / (self.n_h - 1)), self.t0) for i in range(self.n_h) for j in range(self.n_b)]
        self.set_boundary_condition(nodes)
        return nodes


    def create_elements(self):
        return [Element([i * self.n_b + j, i * self.n_b + j + 1, (i + 1) * self.n_b + j + 1, (i + 1) * self.n_b + j]) for i in range(self.n_h - 1) for j in range(self.n_b - 1)]


    def calculate_jacobians(self):
        for element in self.elements:
            element.jacobians = J_matrix_all([self.nodes[id] for id in element.ids])


    def calculate_H_matrices(self):
        for element in self.elements:
            element.H_matrix = H_matrix(element.jacobians)


    def is_boundary_condition(self, element, wall_id):
        if element not in self.elements:
            raise ValueError

        start = wall_id
        if wall_id == len(element.ids) - 1:
            end = 0
        else:
            end = start + 1

        return self.nodes[element.ids[start]].bc and self.nodes[element.ids[end]].bc


    def apply_boundary_conditions(self):
        for element in self.elements: 
            H_BC_matrices = np.zeros((4, 4, 4))
            P_vectors = np.zeros((4, 4, 1))
            for wall_id in range(len(element)):
                if self.is_boundary_condition(element, wall_id):
                    H_BC_matrices[wall_id] = H_matrix_BC(element.jacobians[0], wall_id)
                    P_vectors[wall_id] = P_vector(element.jacobians[0], wall_id)
            element.H_BC_matrices = H_BC_matrices
            element.H_BC_matrix = H_BC_matrices.sum(axis=0)
            element.P_vectors = P_vectors
            element.P_vector = P_vectors.sum(axis=0)


    def update_H_matrices(self):
        for element in self.elements:
            element.H_matrix += element.H_BC_matrix


    def aggregate_H(self):
        aggregated = np.zeros((self.n_n, self.n_n))
        for element in self.elements:
            local_i = list(itertools.product([0, 1, 2, 3], [0, 1, 2, 3]))
            global_i = list(itertools.product(element.ids, element.ids))
            for (l_i, l_j), (g_i, g_j) in zip(local_i, global_i):
                aggregated[g_i, g_j] += element.H_matrix[l_i, l_j]
        return aggregated


    def aggregate_P(self):
        aggregated = np.zeros((self.n_n, 1))
        for element in self.elements:
            for l_i, g_i in zip(range(4), element.ids):
                aggregated[g_i] += element.P_vector[l_i]
        return aggregated


    def aggregate_C(self):
        aggregated = np.zeros((self.n_n, self.n_n))
        for element in self.elements:
            local_i = list(itertools.product([0, 1, 2, 3], [0, 1, 2, 3]))
            global_i = list(itertools.product(element.ids, element.ids))
            for (l_i, l_j), (g_i, g_j) in zip(local_i, global_i):
                aggregated[g_i, g_j] += element.C_matrix[l_i, l_j]
        return aggregated


    def calculate_C_matrices(self):
        for element in self.elements:
            element.C_matrix = C_matrix(element.jacobians[0])

   
    def format_elements(self):
        return '\n'.join([str(element) for element in self.elements])
    

    def __repr__(self):
        return 30 * "=" + f"""\nPrinting info about this grid\n\nNodes:\n{str(self.nodes)}\n\nElements:\n{self.format_elements()}\n\n \
        H global:\n{self.H_global.round(3)}\n\nP global:\n{self.P_global.round(3)}\n\nC global:\n{self.C_global.round(3)}\n\n""" + 30 * "=" 




