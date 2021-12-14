from grid import Grid
from globals import *
import sys


def print_delimiter():
    print("=" * 30)


def main():
    try:
        schema = int(sys.argv[1])
        if schema in (2, 3):
            SCHEMA_N = schema - 2
            print(f"Setting intergation schema... SCHEMA = {SCHEMA_N + 2}")
        else:
            print("Invalid schema. Only schemas for 2 and 3 integration points are supported for now. Using default schema (N = 2)")
    except IndexError:
        print("Integration schema not specified. Using default schema (N = 2).")

    grid = Grid(H, B, N_H, N_B, T_START)    
    for t0, min, max, i in grid.calculate_temperatures(TM_STEP, N_STEPS):
        print(f"After {i * TM_STEP}s: Tmin={min: .3f} Tmax={max: .3f}")


if __name__ == '__main__':
    main()