from jacob import element4_2D


def jacob_all(nodes):
    jacobians = list()

    for xp, yp in zip(element4_2D.get('ksi2w'), element4_2D.get('mu2w')):
        x = y = 0
        for node, nx, ny in zip(nodes, xp, yp):
            x += node.x * nx
            y += node.y * ny
        jacobians.append([[round(x, 7), 0.0], [0.0, round(y, 7)]])

    return jacobians


def jacob(nodes):
    x = y = 0
    for node, nx, ny in zip(nodes, element4_2D.get('ksi2w')[0], element4_2D.get('mu2w')[0]):
        x += node.x * nx
        y += node.y * ny
    return [[[round(x, 6), 0.0], [0.0, round(y, 6)]]] * 4


class Node:

    """
        Represents single point in 2D space.
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y


    def __str__(self):
        return f"[{self.x: .3f}, {self.y: .3f}]" 


    def __repr__(self):
        return self.__str__()
        #return f"Node({self.x}, {self.y})"

class Element:

    """
        Represents single element (rectangle) in 2D space.
    """

    def __init__(self, ids):
        if not isinstance(ids, (list, tuple)):
            raise TypeError
        else:
            self.ids = ids 


    def __str__(self):
        return str(self.ids)


    def __repr__(self):
        return self.__str__()


class Grid:

    """
        Represents single grid that consists of elements and nodes.
    """

    def __init__(self, height, breadth, num_nodes_height, num_nodes_breadth):
        if not isinstance(height, (int, float)) or not isinstance(breadth, (int, float)):
            raise TypeError("Height and breadth have to be numeric values!")
        elif not isinstance(num_nodes_height, int) or not isinstance(num_nodes_breadth, int):
            raise TypeError("Number of nodes along height and breadth has to be integer!")
        else:
            self.height = height
            self.breadth = breadth
            self.num_nodes_breadth = num_nodes_breadth
            self.num_nodes_height = num_nodes_height
            self.num_nodes = self.num_nodes_height * self.num_nodes_breadth
            self.num_elements = (self.num_nodes_height - 1) * (self.num_nodes_breadth - 1)
            self.nodes = self.create_nodes()
            self.elements = self.create_elements()
            self.jacobians = self.create_jacobians()


    def create_nodes(self):
        return [Node(j * (self.breadth / (self.num_nodes_breadth - 1)), i * (self.height / (self.num_nodes_height - 1))) 
                    for i in range(self.num_nodes_height)
                    for j in range(self.num_nodes_breadth)]


    def create_elements(self):
        return [Element([i * self.num_nodes_breadth + j, i * self.num_nodes_breadth + j + 1, (i + 1) * self.num_nodes_breadth + j + 1, (i + 1) * self.num_nodes_breadth + j])
                    for i in range(self.num_nodes_height - 1)
                    for j in range(self.num_nodes_breadth - 1)]


    def create_jacobians(self):
        return [jacob_all([self.nodes[id] for id in element.ids]) for element in self.elements]


    def format_jacobians(self):
        return '\n' + '\n'.join([str(matrix_el) for matrix_el in self.jacobians])
            

    def __str__(self):
        return 30 * "=" + f"""\nPrinting info about this grid\n\nNodes:\n{str(self.nodes)}\n\nElements:\n{str(self.elements)}\n\nJacobians:\n{self.jacobians}\n"""


    def __repr__(self):
        return self.__str__()
        #return f"Grid({self.height},{self.breadth}, {self.num_nodes_height}, {self.num_nodes_breadth}"

