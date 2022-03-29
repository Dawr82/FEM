import sys
from pathlib import Path

from parse import GridDataParser
from grid import Grid
import param


def main_from_data():
    grid = Grid.from_data(
        param.H, param.B, param.N_H, 
        param.N_B, param.T_START)
    for min, max, i in grid.calculate_temperatures(param.TM_STEP, param.N_STEPS):
        print(f"After {i * param.TM_STEP}s: Tmin={min: .3f} Tmax={max: .3f}")


def main_from_file():
    try:
        schema = int(sys.argv[1])
        filename = Path(sys.argv[2])
        if not filename.is_file():
            print(f"This file doesn't exist!")
            sys.exit(-1)
        if schema in (2, 3):
            param.SCHEMA_N = schema - 2
            print(f"Setting intergation schema... SCHEMA_N = {schema}")
        else:
            print("Invalid schema. Only schemas for 2 and 3 \
            integration points are supported for now. \
            Using default schema (N = 2)")
    except IndexError:
        print("Invalid argument list (schema identifier (0 or 1) and filename are required)")
        sys.exit(-1)

    parser = GridDataParser(filename)
    params = parser.params

    param.TM_STEP = params.get("SimulationStepTime", params.get("dt", param.TM_STEP))
    param.N_STEPS = int(params.get("SimulationTime", 1_000_000_000) // param.TM_STEP)
    param.K = params.get("Conductivity", params.get("K", param.K))
    param.ALFA = params.get("Alfa", param.ALFA)
    param.T_AMBIENT = params.get("Tot", param.T_AMBIENT)
    param.T_START = params.get("InitialTemp", param.T_START)
    param.RO = params.get("Density", param.RO)
    param.C = params.get("SpecificHeat", param.K)
   
    grid = Grid.from_file(parser.nodes, parser.elements, parser.bc, param.T_START)

    output = ""

    try:
        for min, max, i in grid.calculate_temperatures(param.TM_STEP, param.N_STEPS):
            curr_output = f"After {i * param.TM_STEP}s: Tmin={min: .3f} Tmax={max: .3f}\n"
            output += curr_output
            print(curr_output, end="")
    except KeyboardInterrupt:
        print("Computation ended.")

    choice = input("\nSave data to a file? Y/N ")

    if choice.upper() == "Y":
        with open(Path(filename.parent) / ("output_" + filename.name) , "w") as file:
            file.write(output)


if __name__ == '__main__':
    main_from_file()
    #main_from_data()