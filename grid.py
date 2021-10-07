class Node:

    """
        Represents single point in 2D space.
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __repr__(self):
        return f"[{self.x}, {self.y}]" 

class Element:

    """
        Represents single element (rectangle) in 2D space.
    """

    def __init__(self, ids):
        if not isinstance(ids, (list, tuple)):
            raise TypeError
        else:
            self.ids = ids   
    def __repr__(self):
        return str(self.ids)

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
    def create_nodes(self):
        return [Node(j * (self.breadth / (self.num_nodes_breadth - 1)), i * (self.height / (self.num_nodes_height - 1))) 
                    for i in range(self.num_nodes_height)
                    for j in range(self.num_nodes_breadth)]
    def create_elements(self):
        return [Element([i * self.num_nodes_breadth + j, i * self.num_nodes_breadth + j + 1, (i + 1) * self.num_nodes_breadth + j, (i + 1) * self.num_nodes_breadth + j + 1])
                    for i in range(self.num_nodes_height - 1)
                    for j in range(self.num_nodes_breadth - 1)]
    def __repr__(self):
        return 30 * "=" + f"""\nPrinting info about this grid\n\nNodes:\n{str(self.nodes)}\n\nElements:\n{str(self.elements)}\n"""

