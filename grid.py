from jacob import element4_2D
import numpy as np

K = 30

def jacob_all(nodes):
    jacobians = list()

    for xp, yp in zip(element4_2D.get('ksi2w'), element4_2D.get('mu2w')):
        x_ksi = x_mu = y_ksi = y_mu = 0
        for node, n_ksi, n_mu in zip(nodes, xp, yp):
            x_ksi += node.x * n_ksi
            y_ksi += node.y * n_ksi
            x_mu += node.x * n_mu
            y_mu += node.y * n_mu
        jacobians.append([[round(x_ksi, 7), round(y_ksi, 7)], [round(x_mu, 7), round(y_mu, 7)]])

    return jacobians


def jacob(nodes):
    x_ksi = x_mu = y_ksi = y_mu = 0
    for node, n_ksi, n_mu in zip(nodes, element4_2D.get('ksi2w')[0], element4_2D.get('mu2w')[0]):
        x_ksi += node.x * n_ksi
        y_ksi += node.y * n_ksi
        x_mu += node.x * n_mu
        y_mu += node.y * n_mu
    #return np.array([x_ksi, y_ksi, x_mu, y_mu]).round(7)
    return [[[round(x_ksi, 6), round(y_ksi, 6)], [round(x_mu, 6), round(y_mu, 6)]]] * 4


def H_matrix(J):
    H_list = list()
    for nx, ny, j in zip(element4_2D.get('ksi2w'), element4_2D.get('mu2w'), J):
        rj = np.linalg.inv(np.array(j))
        dNdx = rj[0, 0] * np.array([nx]) + rj[0, 1] * np.array([nx])
        dNdy = rj[1, 0] * np.array([ny]) + rj[1, 1] * np.array([ny])
        H_list.append((K * (dNdx.T @ dNdx + dNdy.T @ dNdy) * np.linalg.det(j)).round(3).tolist())
    return H_list

   
class Node:

    """
        Represents a single point in 2D space.
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y


    def __repr__(self):
        return f"[{self.x: .3f}, {self.y: .3f}]" 


class Element:

    """
        Represents a single element (rectangle) in 2D space.
    """

    def __init__(self, ids):
        self.ids = ids 
        self.H_matrices = None
        self.jacobians = None
   

    def format_H_matrices(self):
        return '\t'.join([str(matrix) for matrix in self.H_matrices])


    def format_jacobians(self):
        return '\t'.join([str(matrix_el) for matrix_el in self.jacobians])


    def __repr__(self):
        return f"{self.ids}\nJacobians:\n{self.format_jacobians()}\nH matrices:\n{self.format_H_matrices()}\n"


class Grid:

    """
        Represents single grid that consists of elements and nodes.
    """

    def __init__(self, height, breadth, num_nodes_height, num_nodes_breadth):
        self.height = height
        self.breadth = breadth
        self.num_nodes_breadth = num_nodes_breadth
        self.num_nodes_height = num_nodes_height
        self.num_nodes = self.num_nodes_height * self.num_nodes_breadth
        self.num_elements = (self.num_nodes_height - 1) * (self.num_nodes_breadth - 1)
        self.nodes = self.create_nodes()
        self.elements = self.create_elements()
        
        self.calculate_jacobians()
        self.calculate_H_matrices()


    def create_nodes(self):
        return [Node(j * (self.breadth / (self.num_nodes_breadth - 1)), i * (self.height / (self.num_nodes_height - 1))) 
                    for i in range(self.num_nodes_height)
                    for j in range(self.num_nodes_breadth)]


    def create_elements(self):
        return [Element([i * self.num_nodes_breadth + j, i * self.num_nodes_breadth + j + 1, (i + 1) * self.num_nodes_breadth + j + 1, (i + 1) * self.num_nodes_breadth + j])
                    for i in range(self.num_nodes_height - 1)
                    for j in range(self.num_nodes_breadth - 1)]


    def calculate_jacobians(self):
        for element in self.elements:
            element.jacobians = jacob_all([self.nodes[id] for id in element.ids])


    def calculate_H_matrices(self):
        for element in self.elements:
            element.H_matrices = H_matrix(element.jacobians)
    

    def __str__(self):
        return 30 * "=" + f"""\nPrinting info about this grid\n\nNodes:\n{str(self.nodes)}\n\nElements:\n{str(self.elements)}"""


    def __repr__(self):
        return self.__str__()
        #return f"Grid({self.height},{self.breadth}, {self.num_nodes_height}, {self.num_nodes_breadth}"

