from jacob import J_matrix, J_matrix_all, H_matrix, H_matrix_BC, P_vector, C_matrix
from globals import SCHEMA_N
import numpy as np
import itertools

np.set_printoptions(linewidth=200)


class Node:

    def __init__(self, x, y, t_start, bc=0):
        self.x = x
        self.y = y
        self.bc = bc
        self.t_start = t_start

    def __repr__(self):
        return f"Node data:\n\n[{self.x: .3f}, {self.y: .3f}, BC = {self.bc}, T = {self.t_start: .3f}]\n" 


class Element:

    def __init__(self, ids, jacobians=None, H_matrix=None,  H_BC_matrix=None, P_vector=None, C_matrix=None):
        self.ids = ids 
        
        self.jacobians = jacobians
        self.H_matrix = H_matrix
        self.H_BC_matrix = H_BC_matrix
        self.P_vector = P_vector
        self.C_matrix = C_matrix


    def __repr__(self):
        return f"Element data:\n\n{self.ids}\n\nJ matrix:\n\n{self.jacobians}\n\nH matrix:\n\n{self.H_matrix}\n\nH_BC_matrix\n\n{self.H_BC_matrix}\n\nP_vector\n\n{self.P_vector}\n\nC_matrix\n\n{self.C_matrix}\n"


    def __len__(self):
        return len(self.ids)

  
class Grid:

    def __init__(self, height, breadth, num_nodes_height, num_nodes_breadth, t_start):
        self.h = height
        self.b = breadth
        self.n_b = num_nodes_breadth
        self.n_h = num_nodes_height
        self.l_b = self.b / (self.n_b - 1)   # Wysokosc elementu
        self.l_h = self.h / (self.n_h - 1)   # Szerokosc elementu
        self.n_n = self.n_h * self.n_b
        self.n_e = (self.n_h - 1) * (self.n_b - 1)
        self.t_start = t_start

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
        nodes = [Node(j * (self.b / (self.n_b - 1)), i * (self.h / (self.n_h - 1)), self.t_start) for i in range(self.n_h) for j in range(self.n_b)]
        self.set_boundary_condition(nodes)
        return nodes


    def create_elements(self):
        return [Element([i * self.n_b + j, i * self.n_b + j + 1, (i + 1) * self.n_b + j + 1, (i + 1) * self.n_b + j]) for i in range(self.n_h - 1) for j in range(self.n_b - 1)]


    def calculate_jacobians(self):
        for element in self.elements:
            element.jacobians = J_matrix_all([self.nodes[id] for id in element.ids], SCHEMA_N)


    def calculate_H_matrices(self):
        for element in self.elements:
            element.H_matrix = H_matrix(element.jacobians[0], SCHEMA_N)


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
                    if wall_id in(0, 2):
                        H_BC_matrices[wall_id] = H_matrix_BC(self.l_b, wall_id, SCHEMA_N)
                        P_vectors[wall_id] = P_vector(self.l_b, wall_id, SCHEMA_N)
                    elif wall_id in (1, 3):
                        H_BC_matrices[wall_id] = H_matrix_BC(self.l_h, wall_id, SCHEMA_N)
                        P_vectors[wall_id] = P_vector(self.l_h, wall_id, SCHEMA_N)
            element.H_BC_matrix = H_BC_matrices.sum(axis=0)
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
            element.C_matrix = C_matrix(element.jacobians[0], SCHEMA_N)

   
    def format_elements(self):
        return '\n'.join([str(element) for element in self.elements])

    
    def calculate_temperatures(self, dT, n_steps):
        H = self.H_global
        P = self.P_global
        C = self.C_global
        T0 = np.full(P.shape, self.t_start)

        X = H + C/dT
        Y = C/dT

        for i in range(n_steps):
            T0 = np.linalg.inv(X) @ (Y @ T0 + P)
            yield (T0, T0.min(), T0.max(), i + 1)


    def __repr__(self):
        return 30 * "=" + f"""\nPrinting info about this grid\n\nNodes:\n{str(self.nodes)}\n\nElements:\n{self.format_elements()}\n\n \
        H global:\n{self.H_global.round(3)}\n\nP global:\n{self.P_global.round(3)}\n\nC global:\n{self.C_global.round(3)}\n\n""" + 30 * "=" 




