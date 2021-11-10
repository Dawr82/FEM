from jacob import element4_2D
import numpy as np

K = 30

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

   
class Node:

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"[{self.x: .3f}, {self.y: .3f}]" 


class Element:

    def __init__(self, ids, jacobians=None, H_matrices=None):
        self.ids = ids 
        self.jacobians = jacobians
        self.H_matrix = H_matrices


    def __repr__(self):
        return 30 * "=" + f"\n\n{self.ids}\n\nJ matrix:\n\n{self.jacobians}\n\nH matrix:\n\n{self.H_matrix}\n\n"

  
class Grid:

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
            element.jacobians = J_matrix_all([self.nodes[id] for id in element.ids])


    def calculate_H_matrices(self):
        for element in self.elements:
            element.H_matrix = H_matrix(element.jacobians)

    
    def format_elements(self):
        return '\n'.join([str(element) for element in self.elements])
    

    def __repr__(self):
        return 30 * "=" + f"""\nPrinting info about this grid\n\nNodes:\n{str(self.nodes)}\n\nElements:\n{self.format_elements()}\n""" + 30 * "="




def test_H_matrix():
    nodes = [
        Node(0, 0),
        Node(0.025, 0),
        Node(0.025, 0.025),
        Node(0, 0.025)
    ]

    J = J_matrix_all(nodes)
    H = H_matrix(J)
    print(H)

test_H_matrix()