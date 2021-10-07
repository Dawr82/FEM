from grid import *

if __name__ == '__main__':

    # Testing Classes from grid.py here

    # Node
    node = Node(10, 20)
    print(node)

    # Element
    element1 = Element((3,4,1,2))
    element2 = Element([4,5,6,7])

    try:
        element3 = Element({'1' : 3, '2' : 4})
    except TypeError:
        print("Not list-like iterable passed")
    
    try: 
        element4 = Element(5)
    except TypeError:
        print("Single ID passed")

    # Grid
    grid = Grid(0.2, 0.1, 5, 4)
    print(grid.create_nodes())
    print(grid.create_elements())
    print(grid)