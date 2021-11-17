from jacob import element4_2D
from integration import N_SCHEMA_BC_2W
import numpy as np

K = 30
ALFA = 25

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


class Node:

    def __init__(self, x, y, bc=0):
        self.x = x
        self.y = y
        self.bc = bc

    def __repr__(self):
        return f"[{self.x: .3f}, {self.y: .3f}, BC = {self.bc}]" 


class Element:

    def __init__(self, ids, jacobians=None, H_matrix=None, H_BC_matrices=None, H_BC_matrix=None):
        self.ids = ids 
        self.jacobians = jacobians
        self.H_matrix = H_matrix
        self.H_BC_matrices = H_BC_matrices
        self.H_BC_matrix = H_BC_matrix


    def __repr__(self):
        return 30 * "=" + f"\n\n{self.ids}\n\nJ matrix:\n\n{self.jacobians}\n\nH matrix:\n\n{self.H_matrix}\n\n \
H_BC matrices\n\n{self.H_BC_matrices}\n\nH_BC_matrix\n\n{self.H_BC_matrix}\n\n"


    def __len__(self):
        return len(self.ids)

  
class Grid:

    def __init__(self, height, breadth, num_nodes_height, num_nodes_breadth):
        self.h = height
        self.b = breadth
        self.n_b = num_nodes_breadth
        self.n_h = num_nodes_height
        self.n_n = self.n_h * self.n_b
        self.n_e = (self.n_h - 1) * (self.n_b - 1)
        self.nodes = self.create_nodes()
        self.elements = self.create_elements()
  
        self.calculate_jacobians()
        self.calculate_H_matrices()
        self.apply_boundary_conditions()
        #self.update_H_matrices()


    def set_boundary_condition(self, nodes):
        for node in nodes:
            if any((not node.x, not node.y, node.x == self.b, node.y == self.h)):
                node.bc = 1


    def create_nodes(self):
        nodes = [Node(j * (self.b / (self.n_b - 1)), i * (self.h / (self.n_h - 1))) for i in range(self.n_h) for j in range(self.n_b)]
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

        return self.nodes[element.ids[start]].bc and self.nodes[element.ids[end] - 1].bc


    def apply_boundary_conditions(self):
        for element in self.elements: 
            H_BC_matrices = np.zeros((4, 4, 4))
            for wall_id in range(len(element)):
                if self.is_boundary_condition(element, wall_id):
                    H_BC_matrices[wall_id] = H_matrix_BC(element.jacobians[0], wall_id)
            element.H_BC_matrices = H_BC_matrices
            element.H_BC_matrix = H_BC_matrices.sum(axis=0)



    def update_H_matrices(self):
        for element in self.elements:
            element.H_matrix += element.H_BC_matrix
           
   
    def format_elements(self):
        return '\n'.join([str(element) for element in self.elements])
    

    def __repr__(self):
        return 30 * "=" + f"""\nPrinting info about this grid\n\nNodes:\n{str(self.nodes)}\n\nElements:\n{self.format_elements()}\n""" + 30 * "="



