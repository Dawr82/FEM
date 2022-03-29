H = 0.1           # m
B = 0.1           # m
N_H = 4           # nodes per H
N_B = 4           # nodes per B
C = 700           # J/kg*C  - specific heat
K = 25            # W/m*C   - conductivity
RO = 7800         # kg/m3   - density
ALFA = 300        # W/m2K   
T_AMBIENT = 1200  # Temperature of the environment (ambient temperature)
T_START = 100     # Temperature at the beginning of the simulation (starting point for the simulation)
TM_STEP = 50      # Time step (50 seconds)
N_STEPS = 10      # Number of steps
SCHEMA_N = 1      # Indicates integrations schema (Gaussian integration)

# So: We are simulating the process that would take TM_STEP * N_STEPS time in real world.
# We are calculating the temperature in TM_STEP intervals, N_STEPS times.

# Note that this parameters are used only if there is no Solidworks file as an input to the program.
# If there is one, there should be these parameters specified there.

# User of this program is responsible for passing valid SCHEMA_N identifier as an argument to the program.
# If he doesn't do that, defauly schema is applied (1 in this case).