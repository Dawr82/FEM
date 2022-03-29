class GridDataParser:

    def __init__(self, filename, parse=True):
        self.filename = filename
        if parse:
            self.parse()

    
    def parse(self):
        with open(self.filename, "r") as file:
            # parameters
            params = dict()
            while True:
                param_data = file.readline()
                if "*Node" in param_data:
                    break
                data = param_data.rstrip("\n").rsplit(" ", 1)
                params[data[0]] = float(data[1])

            self._params = params

            # nodes
            nodes = list()
            while True:
                node_data = file.readline()
                if "*Element" in node_data:
                    break
                nodes.append([float(n) for n in node_data.lstrip().rstrip("\n").split(", ")[1:]])

            self._nodes = nodes

            # elements
            elements = list()            
            while True:
                element_data = file.readline()
                if "*BC" in element_data:
                    break
                elements.append([int(n) - 1 for n in element_data.lstrip().rstrip("\n").split(", ")[1: ]])

            self._elements = elements

            # bc
            self._bc = [int(n) - 1 for n in file.readline().rstrip("\n").split(", ")]


    @property
    def params(self):
        return self._params

    
    @property
    def nodes(self):
        return self._nodes


    @property
    def elements(self):
        return self._elements

    @property
    def bc(self):
        return self._bc