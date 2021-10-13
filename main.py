from grid import *
from integration import *
import pprint

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

    print("\n\nIntegrals:\n")

    f1 = lambda x : 5*x**2 + 3*x + 6
    f1_1N = gaussian_integration1D(f1, 1)
    f1_2N = gaussian_integration1D(f1, 2) 
    print(f1_1N, f1_2N)

    f2 = lambda x, y : 5 * x**2 * y**2 + 3*x*y + 6
    f2_1N = gaussian_integration2D(f2, 1)
    f2_2N = gaussian_integration2D(f2, 2)
    print(f2_1N, f2_2N)