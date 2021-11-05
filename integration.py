from math import sqrt
from inspect import signature

class IntegrationError(Exception):
    pass

SCHEMA = [
    [[1, 1],[-1/sqrt(3), 1/sqrt(3)]],                                                            # 2 w
    [[5/9, 8/9, 5/9],[-sqrt(3/5), 0, sqrt(3/5)]],                                                # 3 w
    [[0.347855, 0.652145, 0.652145,  0.347855],[-0.861136, -0.339981, 0.339981, 0.861136]],      # 4 w
]

def integral_gauss(f, N):
    n_vars = len(signature(f).parameters)
    result = 0

    try:
        if n_vars == 1:
            for i in range(N + 2):
                result += f(SCHEMA[N][1][i]) * SCHEMA[N][0][i]
        elif n_vars == 2:
            for i in range(N + 2):
                for j in range(N + 2):
                    result += f(SCHEMA[N][1][i], SCHEMA[N][1][j]) * SCHEMA[N][0][i] * SCHEMA[N][0][j]
        else:
            raise IntegrationError("Function has to have 1 or 2 variables.")
    except IndexError:
        print("SCHEMA doesn't support that N.")

    return result



