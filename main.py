from parse import GridDataParser
from grid import Grid
import constants

import sys
from pathlib import Path


def main1():
    grid = Grid.from_data(constants.H, constants.B, constants.N_H, constants.N_B, constants.T_START)
    for min, max, i in grid.calculate_temperatures(constants.TM_STEP, constants.N_STEPS):
        print(f"After {i * constants.TM_STEP}s: Tmin={min: .3f} Tmax={max: .3f}")


def main():
    try:
        schema = int(sys.argv[1])
        filename = Path(sys.argv[2])
        if not filename.is_file():
            print(f"This file doesn't exist!")
            sys.exit()
        if schema in (2, 3):
            SCHEMA_N = schema - 2
            print(f"Setting intergation schema... SCHEMA = {schema}")
        else:
            print("Invalid schema. Only schemas for 2 and 3 integration points are supported for now. Using default schema (N = 2)")
    except IndexError:
        print("Invalid argument list (schema identifier (0 or 1) and filename are required)")
        sys.exit(0)

    parser = GridDataParser(filename)

    params = parser.params

    TM_STEP = params.get("SimulationStepTime", params.get("dt", constants.TM_STEP))
    N_STEPS = int(params.get("SimulationTime", 1_000_000_000) // TM_STEP)
    K = params.get("Conductivity", params.get("K", constants.K))
    ALFA = params.get("Alfa", constants.ALFA)
    T_AMBIENT = params.get("Tot", constants.T_AMBIENT)
    T_START = params.get("InitialTemp", constants.T_START)
    RO = params.get("Density", constants.RO)
    C = params.get("SpecificHeat", constants.K)
   
    grid = Grid.from_file(parser.nodes, parser.elements, parser.bc, T_START)

    output = ""

    try:
        for min, max, i in grid.calculate_temperatures(TM_STEP, N_STEPS):
            curr_output = f"After {i * TM_STEP}s: Tmin={min: .3f} Tmax={max: .3f}\n"
            output += curr_output
            print(curr_output, end="")
    except KeyboardInterrupt:
        print("Computation ended.")

    choice = input("\nSave data to file? Y/N ")

    if choice.upper() == "Y":
        with open(Path(filename.parent) / ("output_" + filename.name) , "w") as file:
            file.write(output)


if __name__ == '__main__':
    #main()
    main1()